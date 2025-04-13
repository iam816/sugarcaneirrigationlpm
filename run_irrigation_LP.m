function results = run_irrigation_LP(budget, w_min, threshold, Cw, Ce, Cs, min_zones, max_zones)
    % Load the soil moisture data
    data = readtable('september2014_soilmoisture.csv');

    num_days = height(data);

    results = table('Size', [num_days, 6], ...
        'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'string'}, ...
        'VariableNames', {'Date','TotalCost','WaterUsed','EnergyUsed','ZonesIrrigated','ActivatedZones'});

    for day = 1:num_days
        date_str = string(data.Date(day));
        moisture = [data.Zone1(day), data.Zone2(day), data.Zone3(day), data.Zone4(day)];
        deficit = max(1 - moisture, 0.01);
        Cw_dynamic = Cw ./ deficit;

        f = [Cw_dynamic, Ce*ones(1,4), Cs*ones(1,4)];
        intcon = 9:12;
        lb = zeros(12,1);
        ub = [9500 * ones(1,4), 21.11 * ones(1,4), ones(1,4)];

        % Constraints
        A1 = f;                         b1 = budget;
        A2 = [1 1 1 1, zeros(1,8)];     b2 = 38000;
        A3 = [eye(4), zeros(4,4), -9500*eye(4)]; b3 = zeros(4,1);
        A4 = [zeros(4,4), eye(4), -21.11*eye(4)]; b4 = zeros(4,1);
        A5 = [zeros(1,8), ones(1,4)];   b5 = max_zones;

        A6 = [-eye(4), zeros(4,4), w_min * eye(4)];
        b6 = zeros(4,1);

        A = [A1; A2; A3; A4; A5; A6];
        b = [b1; b2; b3; b4; b5; b6];

        % Min number of activated zones: x1 + x2 + ... ≥ min_zones → -x ≤ -min_zones
        A = [A; zeros(1,8), -ones(1,4)];
        b = [b; -min_zones];

        Aeq = [-0.002222 * eye(4), eye(4), zeros(4,4)];
        beq = zeros(4,1);

        options = optimoptions('intlinprog','Display','off');
        [x, fval, exitflag] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, options);

        if exitflag == 1
            w_used = sum(x(1:4));
            e_used = sum(x(5:8));
            activated = find(round(x(9:12)) == 1);
            results(day,:) = {date_str, fval, w_used, e_used, length(activated), join("Q"+activated, ', ')};
        else
            results(day,:) = {date_str, NaN, NaN, NaN, 0, "No solution"};
        end
    end
end
