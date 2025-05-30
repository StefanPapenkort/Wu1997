% -------------------- Implementation Huxley Model ------------------------

% spatial domain: x (-2,2) [mm]
% temporal domain: t (0,0.2) [s]

h = 1;
xi = [-2,2]/h; % normalized spatial coordinate

% introduce spatial discritization
xi_i = linspace(xi(1),xi(2),20);
numPts = length(xi_i);

% determine Gauss Points and corresponding intervals
GaussPoints = zeros((numPts-1)*2,5);
numGaussPts = size(GaussPoints,1);
for j = 1:numPts-1
    gp1 = xi_i(j) + ((3 - sqrt(3))/6)*(xi_i(j+1)-xi_i(j));
    gp2 = xi_i(j) + ((3 + sqrt(3))/6)*(xi_i(j+1)-xi_i(j));
    GaussPoints(2*j-1,:) = [gp1, xi_i(j), xi_i(j+1),j,j+1]; % [GP coord, P1 coord, P2 coord, P1 ID, P2 ID]
    GaussPoints(2*j,:) = [gp2, xi_i(j), xi_i(j+1),j,j+1];   % [GP coord, P1 coord, P2 coord, P1 ID, P2 ID]
end



%% --------------------- R U N   S I M U L A T I O N -----------------------

t_max = 0.3;    % simulation time in [s]
dt = 0.001;       % step size [s]

numSteps = t_max/dt;

% initialize fields that store da_i/dt and db_i/dt for 4 previous time steps
dAdt_prev = zeros(numPts,4);
dBdt_prev = zeros(numPts,4);

% initialize a_i and b_i (here a_prev and b_prev are a0 and b0, respectively)
a_prev = zeros(numPts,1);
Idx = xi_i > 0 & xi_i < 1;
a_prev(Idx) = 43.3/(43.3 + 10);
b_prev = zeros(numPts,1);

x0 = [a_prev; a_prev*0; b_prev; b_prev*0]; % starting value for simulation

% arrays for storing the solution
A_sol = zeros(numSteps+1,numPts); A_sol(1,:) = a_prev';
B_sol = zeros(numSteps+1,numPts); B_sol(1,:) = b_prev';

for step = 1:numSteps

    f = @(x) AdamsMoulton(x, a_prev, dAdt_prev, b_prev, dBdt_prev, GaussPoints, step, dt); % function of dummy variable x
    %OPTIONS = optimoptions("fsolve",'Algorithm','levenberg-marquardt','FiniteDifferenceType','central');
    [x, FVAL, info] = fsolve(f,x0); %,OPTIONS);

    % update variables
    a       = x(1:numPts);
    a_prime = x(numPts+1:2*numPts);
    b       = x(2*numPts+1:3*numPts);
    b_prime = x(3*numPts+1:4*numPts);


    a_prev = a;
    dAdt_prev(:,2:4) = dAdt_prev(:,1:3);
    dAdt_prev(:,1) = a_prime;

    b_prev = b;
    dBdt_prev(:,2:4) = dBdt_prev(:,1:3);
    dBdt_prev(:,1) = b_prime;

    x0 = x;   % update starting value for next time step

    % store results in solution arrays
    A_sol(step+1,:) = a_prev';
    B_sol(step+1,:) = b_prev';
end

%% evaluate and visualize n(x,t)
xi_eval = -2:0.01:2;
xi_eval_mat = zeros(length(xi_eval),5);
for i = 1:length(xi_eval)
    flag = false;
    j = 1;
    while flag == false
        if xi_eval(i) >= xi_i(j) && xi_eval(i) <= xi_i(j+1)
            xi_eval_mat(i,:) = [xi_eval(i), xi_i(j), xi_i(j+1), j, j+1];
            flag = true;
        end
        j=j+1;
    end
end

n_sol = zeros(numSteps+1,size(xi_eval_mat,1)); % create empty array
for i=1:numSteps+1
    for j = 1:size(xi_eval_mat,1)
        n_sol(i,j) = eval_n(A_sol(i,:)',B_sol(i,:)',xi_eval_mat(j,:));
    end
end


t = (0:dt:t_max)';
x = xi_eval*h;


[X,T] = meshgrid(x,t);
figure(1)
surf(T,X,n_sol)
xlabel('t [s]')
ylabel('x [mm]')
zlabel('n(x,t) [-]')
zlim([-1,1])

