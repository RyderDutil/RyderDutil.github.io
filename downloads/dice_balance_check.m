%% Name: Ryder
% Fixed face arrangement:
%   +X = PLA    (1 pip)
%   -X = Brass  (6 pips)
%   +Y = Steel  (2 pips)
%   -Y = Acrylic(5 pips)
%   +Z = Wood   (3 pips)
%   -Z = Al     (4 pips)

clc; clear; close all;

%% Constants
SIDE   = 100;    % mm
PIP_M  = 1.5;    % mm, min clearance between pip edges
EDGE_M = 5.0;    % mm, min clearance from pip edge to face edge
D_MIN  = 5.0;    % mm, min pip diameter
D_MAX_PLA = 60;  % mm, max pip diameter for PLA (single center pip)
D_MAX     = 40;  % mm, max pip diameter for all other faces

%% Densities (g/mm^3)
RHO.PLA     = 0.0013;
RHO.Brass   = 0.0085;
RHO.Steel   = 0.008;
RHO.Acrylic = 0.0012;
RHO.Wood    = 0.00063;
RHO.Al      = 0.0027;

%% Discrete thickness options (mm)
DISC.Brass   = [0.2, 0.8, 2.5, 6.5];
DISC.Steel   = [0.2, 0.9, 1.5, 2.0];
DISC.Acrylic = [1.6, 4.3, 4.4, 5.3, 7.9];
DISC.Wood    = [1.2, 2.4, 3.6, 4.8, 6.0, 7.2];
DISC.Al      = [1.7, 3.0, 4.5, 6.6];

%% Build all thickness combos
combos = build_combos(DISC);
n_combos = size(combos, 1);
fprintf('Total discrete thickness combos: %d\n\n', n_combos);

%% Solver options
opts = optimset('Display','off','MaxIter',2000,'MaxFunEvals',50000,'TolFun',1e-9,'TolX',1e-9);
TOL_CG   = 2.0;
TOL_M    = 3.0;
TOL_PROD = 0.2;

%% Main sweep
solutions = {};
best_all  = {};  % track best regardless of pass/fail

fprintf('Running... (this may take a few minutes)\n');
t_start = tic;

for ci = 1:n_combos
    tBrass   = combos(ci,1);
    tSteel   = combos(ci,2);
    tAcrylic = combos(ci,3);
    tWood    = combos(ci,4);
    tAl      = combos(ci,5);

    % Solve PLA analytically from X-axis COM balance
    tPLA = solve_PLA(tBrass, RHO);
    if isnan(tPLA) || tPLA <= 0 || tPLA > 8
        continue;
    end

    ix = SIDE - tPLA - tBrass;
    iy = SIDE - tSteel - tAcrylic;
    if ix <= 0 || iy <= 0 || (SIDE - tWood - tAl) <= 0
        continue;
    end

    T.tPLA=tPLA; T.tBrass=tBrass; T.tSteel=tSteel;
    T.tAcrylic=tAcrylic; T.tWood=tWood; T.tAl=tAl;

    % Initial params: [RC1,d1, RC2,d2, ..., RC6,d6]
    P0 = initial_params(T, SIDE, EDGE_M, PIP_M, D_MIN, D_MAX, D_MAX_PLA);

    % Multi-start: 6 starts exploring different RC/d combinations
    starts = {
        P0;
        P0 .* repmat([0.4, 0.9], 1, 6);
        P0 .* repmat([0.9, 0.2], 1, 6);
        P0 .* repmat([0.5, 0.5], 1, 6);
        P0 .* repmat([0.7, 0.8], 1, 6);
        P0 .* repmat([0.8, 0.4], 1, 6);
    };

    best_f = Inf;
    best_x = P0;
    for si = 1:length(starts)
        p0 = starts{si};
        [x, fval] = fminsearch(@(P) objective(P, T, SIDE, EDGE_M, PIP_M, D_MIN, D_MAX, D_MAX_PLA, RHO, TOL_CG, TOL_M, TOL_PROD), p0, opts);
        if fval < best_f
            best_f = fval;
            best_x = x;
        end
    end

    % Evaluate final result
    r = evaluate(best_x, T, SIDE, EDGE_M, PIP_M, D_MIN, D_MAX, D_MAX_PLA, RHO);

    entry.T = T;
    entry.r = r;

    % Track best regardless
    if isempty(best_all) || r.score < best_all{1}.r.score
        best_all = [{entry}, best_all];
    end

    if r.cg_err <= TOL_CG && r.m_err <= TOL_M && r.prod_err <= TOL_PROD
        solutions{end+1} = entry; %#ok<SAGROW>
        fprintf('  SOLUTION #%d: dI=%.3f%% CG=%.3f%% Prod=%.3f%%  PLA=%.3f Br=%.1f St=%.1f Ac=%.1f Wd=%.1f Al=%.1f\n', ...
            length(solutions), r.m_err, r.cg_err, r.prod_err, ...
            tPLA, tBrass, tSteel, tAcrylic, tWood, tAl);
    end

    if mod(ci, 50) == 0
        fprintf('  Progress: %d/%d combos  solutions found: %d\n', ci, n_combos, length(solutions));
    end
end

elapsed = toc(t_start);
fprintf('\nDone. %.1f seconds. %d/%d combos. %d solutions.\n\n', elapsed, n_combos, n_combos, length(solutions));

%% Sort and display results
if ~isempty(solutions)
    % Sort by moment error + cg error
    scores = cellfun(@(s) s.r.m_err + s.r.cg_err + s.r.prod_err*10, solutions);
    [~, idx] = sort(scores);
    solutions = solutions(idx);

    fprintf('=== TOP SOLUTIONS ===\n');
    n_show = min(10, length(solutions));
    for i = 1:n_show
        print_solution(solutions{i}, i, SIDE);
    end
else
    fprintf('No solutions within all tolerances. Best found:\n\n');
    if ~isempty(best_all)
        % Sort best_all
        scores = cellfun(@(s) s.r.score, best_all);
        [~, idx] = sort(scores);
        best_all = best_all(idx);
        for i = 1:min(5, length(best_all))
            print_solution(best_all{i}, i, SIDE);
        end
    end
end

%% LOCAL FUNCTIONS

function combos = build_combos(DISC)
    % Cross product of all discrete thickness options
    combos = [];
    for tBr = DISC.Brass
    for tSt = DISC.Steel
    for tAc = DISC.Acrylic
    for tWd = DISC.Wood
    for tAl = DISC.Al
        combos(end+1,:) = [tBr, tSt, tAc, tWd, tAl]; %#ok<AGROW>
    end
    end
    end
    end
    end
end

function tPLA = solve_PLA(tBrass, RHO)
    % X-axis COM balance: rho_PLA * t * (SIDE/2 - t/2) = rho_Brass * tBrass * (SIDE/2 - tBrass/2)
    % (rho_PLA/2)*t^2 - (rho_PLA*SIDE/2)*t + C = 0
    SIDE = 100;
    C = RHO.Brass * tBrass * (SIDE/2 - tBrass/2);
    a = RHO.PLA / 2;
    b = -RHO.PLA * SIDE / 2;
    disc = b^2 - 4*a*C;
    if disc < 0
        tPLA = NaN; return;
    end
    roots = [(-b - sqrt(disc))/(2*a), (-b + sqrt(disc))/(2*a)];
    valid = roots(roots > 0 & roots <= 8);
    if isempty(valid)
        tPLA = NaN;
    else
        tPLA = min(valid);
    end
end

function offsets = pip_offsets(slot_id, RC)
    % Global coordinate pip offsets from face center (= global origin)
    switch slot_id
        case 1  % PLA: 1 pip center
            offsets = [0, 0, 0];
        case 2  % Brass: 6 pips at 60 deg in Y-Z plane
            angles = (0:5)' * pi/3;
            offsets = [zeros(6,1), RC*cos(angles), RC*sin(angles)];
        case 3  % Steel: 2 pips along X
            offsets = [-RC, 0, 0; RC, 0, 0];
        case 4  % Acrylic: 5 pips center + cardinal in X-Z
            offsets = [0,0,0; -RC,0,0; RC,0,0; 0,0,-RC; 0,0,RC];
        case 5  % Wood: 3 pips along X
            offsets = [-RC,0,0; 0,0,0; RC,0,0];
        case 6  % Al: 4 pips cardinal in X-Y
            offsets = [-RC,0,0; RC,0,0; 0,-RC,0; 0,RC,0];
    end
end

function RC_max = max_RC(slot_id, ix, iy, EDGE_M)
    SIDE = 100;
    switch slot_id
        case {1,2}  % X faces: SIDE x SIDE
            RC_max = SIDE/2 - EDGE_M;
        case {3,4}  % Y faces: ix x SIDE
            RC_max = min(ix/2, SIDE/2) - EDGE_M;
        case {5,6}  % Z faces: ix x iy
            RC_max = min(ix/2, iy/2) - EDGE_M;
    end
    RC_max = max(RC_max, 1);
end

function d_max = max_pip_d(slot_id, RC, ix, iy, EDGE_M, PIP_M, D_MAX, D_MAX_PLA)
    SIDE = 100;
    switch slot_id
        case {1,2}; hX = SIDE/2; hY = SIDE/2;
        case {3,4}; hX = ix/2;   hY = SIDE/2;
        case {5,6}; hX = ix/2;   hY = iy/2;
    end

    % Outermost pip reach per face
    switch slot_id
        case 1; oX=0;  oY=0;   % PLA center
        case 2; oX=RC; oY=RC;  % Brass 60 deg ring
        case 3; oX=RC; oY=0;   % Steel along X
        case 4; oX=RC; oY=RC;  % Acrylic cardinal X and Z
        case 5; oX=RC; oY=0;   % Wood along X
        case 6; oX=RC; oY=RC;  % Al cardinal X and Y
    end

    if slot_id == 1
        edge_clear = min(hX, hY) - EDGE_M;
    else
        edge_clear = min(hX - oX - EDGE_M, hY - oY - EDGE_M);
    end

    % Min pip-to-pip spacing
    switch slot_id
        case 1; spacing = Inf;
        case 2; spacing = RC;            % 60 deg: chord = RC
        case 3; spacing = 2*RC;          % 2 pips at +-RC
        case 4; spacing = RC;            % center to outer = RC
        case 5; spacing = RC;            % adjacent in column = RC
        case 6; spacing = RC*sqrt(2);    % cardinal: adjacent dist = RC*sqrt(2)
    end

    if isinf(spacing)
        if slot_id == 1
            spacing_clear = D_MAX_PLA;
        else
            spacing_clear = D_MAX;
        end
    else
        spacing_clear = spacing - PIP_M;
    end

    abs_cap = D_MAX;
    if slot_id == 1; abs_cap = D_MAX_PLA; end

    d_max = max(0.5, min([abs_cap, 2*edge_clear, spacing_clear]));
end

function P0 = initial_params(T, SIDE, EDGE_M, PIP_M, D_MIN, D_MAX, D_MAX_PLA)
    ix = SIDE - T.tPLA - T.tBrass;
    iy = SIDE - T.tSteel - T.tAcrylic;
    P0 = zeros(1, 12);
    for i = 1:6
        RC_max = max_RC(i, ix, iy, EDGE_M);
        RC     = max(1, RC_max * 0.5);
        d_max  = max_pip_d(i, RC, ix, iy, EDGE_M, PIP_M, D_MAX, D_MAX_PLA);
        P0(i*2-1) = RC;
        P0(i*2)   = max(D_MIN, d_max * 0.5);
    end
end

function CP = clamp_params(P, T, SIDE, EDGE_M, PIP_M, D_MIN, D_MAX, D_MAX_PLA)
    ix = SIDE - T.tPLA - T.tBrass;
    iy = SIDE - T.tSteel - T.tAcrylic;
    CP = P;
    for i = 1:6
        RC_max   = max_RC(i, ix, iy, EDGE_M);
        RC       = max(0.5, min(RC_max, CP(i*2-1)));
        CP(i*2-1)= RC;
        d_max    = max_pip_d(i, RC, ix, iy, EDGE_M, PIP_M, D_MAX, D_MAX_PLA);
        CP(i*2)  = max(D_MIN, min(d_max, CP(i*2)));
    end
end

function r = evaluate(P, T, SIDE, EDGE_M, PIP_M, D_MIN, D_MAX, D_MAX_PLA, RHO)
    CP = clamp_params(P, T, SIDE, EDGE_M, PIP_M, D_MIN, D_MAX, D_MAX_PLA);

    ix = SIDE - T.tPLA - T.tBrass;
    iy = SIDE - T.tSteel - T.tAcrylic;

    % Face plate centroid offsets (asymmetric bounding correction)
    cx_ip = (T.tBrass   - T.tPLA)    / 2;
    cy_ip = (T.tAcrylic - T.tSteel)  / 2;

    RHO_arr = [RHO.PLA, RHO.Brass, RHO.Steel, RHO.Acrylic, RHO.Wood, RHO.Al];
    TH      = [T.tPLA, T.tBrass, T.tSteel, T.tAcrylic, T.tWood, T.tAl];
    W_arr   = [SIDE, SIDE, ix,   ix,   ix, ix  ];
    H_arr   = [SIDE, SIDE, SIDE, SIDE, iy, iy  ];
    AXIS    = 'xxyyzz';
    SIGN    = [1,-1,1,-1,1,-1];

    tm=0; mx=0; my=0; mz=0;
    Ixx=0; Iyy=0; Izz=0; Ixy=0; Ixz=0; Iyz=0;

    for i = 1:6
        rho = RHO_arr(i);
        w   = W_arr(i); h = H_arr(i); th = TH(i);
        sm  = rho * w * h * th;
        cn  = SIGN(i) * (SIDE/2 - th/2);

        switch AXIS(i)
            case 'x'; cx=cn;    cy=0;    cz=0;
            case 'y'; cx=cx_ip; cy=cn;   cz=0;
            case 'z'; cx=cx_ip; cy=cy_ip;cz=cn;
        end

        % Solid plate centroidal MOI
        switch AXIS(i)
            case 'x'
                ixx_c=sm*(w^2+h^2)/12; iyy_c=sm*(th^2+h^2)/12; izz_c=sm*(th^2+w^2)/12;
            case 'y'
                ixx_c=sm*(th^2+h^2)/12; iyy_c=sm*(w^2+h^2)/12; izz_c=sm*(th^2+w^2)/12;
            case 'z'
                ixx_c=sm*(th^2+h^2)/12; iyy_c=sm*(th^2+w^2)/12; izz_c=sm*(w^2+h^2)/12;
        end

        tm=tm+sm; mx=mx+sm*cx; my=my+sm*cy; mz=mz+sm*cz;
        Ixx=Ixx+ixx_c+sm*(cy^2+cz^2);
        Iyy=Iyy+iyy_c+sm*(cx^2+cz^2);
        Izz=Izz+izz_c+sm*(cx^2+cy^2);
        Ixy=Ixy+sm*cx*cy;
        Ixz=Ixz+sm*cx*cz;
        Iyz=Iyz+sm*cy*cz;

        RC  = CP(i*2-1);
        d   = CP(i*2);
        rr  = d/2;
        hm  = rho * pi * rr^2 * th;

        switch AXIS(i)
            case 'x'; hixx=hm*rr^2/2; hiyy=hm*(3*rr^2+th^2)/12; hizz=hm*(3*rr^2+th^2)/12;
            case 'y'; hixx=hm*(3*rr^2+th^2)/12; hiyy=hm*rr^2/2; hizz=hm*(3*rr^2+th^2)/12;
            case 'z'; hixx=hm*(3*rr^2+th^2)/12; hiyy=hm*(3*rr^2+th^2)/12; hizz=hm*rr^2/2;
        end

        offsets = pip_offsets(i, RC);
        for j = 1:size(offsets,1)
            hx=cx+offsets(j,1); hy=cy+offsets(j,2); hz=cz+offsets(j,3);
            dIxx=hixx+hm*(hy^2+hz^2);
            dIyy=hiyy+hm*(hx^2+hz^2);
            dIzz=hizz+hm*(hx^2+hy^2);
            tm=tm-hm; mx=mx-hm*hx; my=my-hm*hy; mz=mz-hm*hz;
            Ixx=Ixx-dIxx; Iyy=Iyy-dIyy; Izz=Izz-dIzz;
            Ixy=Ixy-hm*hx*hy;
            Ixz=Ixz-hm*hx*hz;
            Iyz=Iyz-hm*hy*hz;
        end
    end

    gx=mx/tm; gy=my/tm; gz=mz/tm;
    cg_err = sqrt(gx^2+gy^2+gz^2) / SIDE * 100;

    I_min = min([Ixx,Iyy,Izz]);
    I_max = max([Ixx,Iyy,Izz]);
    m_err = (I_max - I_min) / I_min * 100;

    prod_max = max([abs(Ixy),abs(Ixz),abs(Iyz)]);
    prod_err = prod_max / I_min * 100;

    r.cg_err=cg_err; r.m_err=m_err; r.prod_err=prod_err;
    r.Ixx=Ixx; r.Iyy=Iyy; r.Izz=Izz;
    r.Ixy=Ixy; r.Ixz=Ixz; r.Iyz=Iyz;
    r.gx=gx; r.gy=gy; r.gz=gz; r.tm=tm;
    r.CP=CP;
    r.score = cg_err + m_err + prod_err*10;
end

function f = objective(P, T, SIDE, EDGE_M, PIP_M, D_MIN, D_MAX, D_MAX_PLA, RHO, TOL_CG, TOL_M, TOL_PROD)
    r = evaluate(P, T, SIDE, EDGE_M, PIP_M, D_MIN, D_MAX, D_MAX_PLA, RHO);
    cg_pen   = max(0, r.cg_err   - TOL_CG)   * 30;
    m_pen    = max(0, r.m_err    - TOL_M)    * 30;
    prod_pen = max(0, r.prod_err - TOL_PROD) * 30;
    f = r.cg_err + r.m_err + r.prod_err*10 + cg_pen + m_pen + prod_pen;
end

function print_solution(entry, rank, ~)
    T = entry.T; r = entry.r; CP = r.CP;
    FACE_NAMES = {'+X PLA','-X Brass','+Y Steel','-Y Acrylic','+Z Wood','-Z Al'};
    fprintf('\nSolution #%d\n', rank);
    fprintf('  dI=%.4f%%  CG=%.4f%%  Prod=%.4f%%\n', r.m_err, r.cg_err, r.prod_err);
    fprintf('  Ixx=%.2f  Iyy=%.2f  Izz=%.2f  (g*mm^2)\n', r.Ixx, r.Iyy, r.Izz);
    fprintf('  Ixy=%.3f  Ixz=%.3f  Iyz=%.3f\n', r.Ixy, r.Ixz, r.Iyz);
    fprintf('  gx=%.4f  gy=%.4f  gz=%.4f mm  mass=%.3fg\n', r.gx, r.gy, r.gz, r.tm);
    fprintf('  Thicknesses (mm): PLA=%.3f  Brass=%.1f  Steel=%.1f  Acrylic=%.1f  Wood=%.1f  Al=%.1f\n', ...
        T.tPLA, T.tBrass, T.tSteel, T.tAcrylic, T.tWood, T.tAl);
    fprintf('  Hole geometry:\n');
    for i = 1:6
        RC = CP(i*2-1); d = CP(i*2);
        fprintf('    %s: CC_diam=%.3fmm  pip_diam=%.3fmm\n', FACE_NAMES{i}, RC*2, d);
    end
end