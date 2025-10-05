N = 480;
C = 46;
S = 39;
CS = 7;

sim_CS = [];
for sim = 1:10000
    sim_C = randperm(N) <= C;
    sim_S = randperm(N) <= S;
    sim_CS(sim) = sum(sim_C .* sim_S);
end

figure;
hist(sim_CS, 10); 
v = axis; hold all; 
plot(CS*[1 1], v(3:4), 'r-');

p = mean(sim_CS >= CS);
title(sprintf('Distrib. of CS, p=%4.3f', p));

%% color vs. axis

N = 480;
C = 46;
A = 111;
CA = 12;

sim_CA = [];
for sim = 1:10000,
    sim_C = (randperm(N) <= C);
    sim_A = (randperm(N) <= A);
    sim_CA(sim) = sum(sim_C .* sim_A);
end

figure;
hist(sim_CA, 10); 
v = axis; hold all; 
plot(CA*[1 1], v(3:4), 'r-');

p = mean(sim_CA >= CA);
title(sprintf('Distrib. of CA, p=%4.3f', p));

%% color vs. shape

N = 480;
S = 39;
A = 111;
SA = 15;

sim_SA = [];
for sim = 1:10000,
    sim_S = (randperm(N) <= S);
    sim_A = (randperm(N) <= A);
    sim_SA(sim) = sum(sim_S .* sim_A);
end

figure;
hist(sim_SA, 10); 
v = axis; hold all; 
plot(SA*[1 1], v(3:4), 'r-');

p = mean(sim_SA >= SA);
title(sprintf('Distrib. of SA, p=%4.3f', p));