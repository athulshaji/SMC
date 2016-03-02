clear;clc;
addpath('../particles_2_1_Matlab')

%% Initialize a particle system object with parameters:
gravity = [0,0,0];                                  % acceleration due to gravity: [0, 0, 0]
drag = 0;                                               % aerodynamic resistance (drag): 0
limits = 1;                                             % axes limits: +/- 1
Particle_System = particle_system(gravity,drag,limits); %

view(3)                                                 % see the plot in 3D
scale = 1;                                              % scale for the plot
axis(scale*[-1 1 -1 1 -1 1])                            % plot will be in this range
%% particle creation with parameters
mass = 42;
initialposition1 = [0,0,0];
velocity = [0,0,0];
fixed_free = false;
lifespan = inf;
Particle1 = particle (Particle_System, mass, initialposition1, velocity, fixed_free, lifespan);

initialposition2 = initialposition1 + 0.5;
Particle2 = particle (Particle_System, mass, initialposition2, velocity, fixed_free, lifespan);
% particle is automatically included in the Particle system
set(Particle2.graphics_handle , 'color', 'blue'); % change properties after creation

%% Spring creation
% A default spring between two already existing particles, Particle_1 and 
% Particle_2, is created and incorporated into the particle system via:
restlength = .2;
springstrength = 100;
dampingcoefficient = 100;
Spring = spring(Particle_System , Particle1 , Particle2,restlength,springstrength,dampingcoefficient);

%% Attraction Creation with parameters
% - gravitational_attraction_strength (which can be made negative for a repulsion) 
% - minimum distance, below which the attraction force will not grow any further. 
% gravattstrg = 0.01;
% mindist = 0.01;
% Attraction = attraction(Particle_System , Particle1 , Particle2);

%% Simulation
step_time = 0.01;
for i = 1:inf
    Particle_System.advance_time(step_time);
end




