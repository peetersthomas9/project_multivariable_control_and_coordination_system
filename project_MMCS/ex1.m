%==========================================================================
%   TP :            Case study: Exercse 1
%   Contact:        ezequiel.gonzalezdebada@epfl.ch
%==========================================================================
classdef ex1
    %Class gathering the solutions of exercise 1. 
    methods (Static)
        %
        function varargout = getSystemParameters
            % PARAMETERS = getSystemParameters() returns a 5-elements
            % column vector containing the value of the system parameters 
            % and the linearization point. Specifically it should contain 
            % (in the presented order):
            %   - k : value of curvature of the reference path.
            %   - car_length [m]: car's length.  
            %   - sigma_v : coefficient characterizing the dynamic of the
            %   actuator tracking the speed reference. 
            %   - sigma_S : coefficient characterizing the dynamic of the
            %   actuator tracking the steering wheel's reference position. 
            %   - spReg : value of the speed the vehicle should drive at. 

            varargout = {[1*10^(-10),4,1,5,5]};               
        end
        %
        function varargout = getLinealModelArrays(parameters)
            % [A,B,C,D] = getLinealModelArrays(PARAMETERS) returns the
            % matrices A,B,C,D characterizing the continuous-time linear
            % model of the system. The input *parameters* corresponds to
            % the output of the method *getSystemParameters*.
            
            %%- Calculate arrays A,B,C,D
            
            % Method to get A,B,C,D-----------------------------
%             syms x1 x2 x3 x4 x5 % state variable
%             syms u1 u2 %input
%             syms t %time vector
%             
%             % state space equation :
%             x1d=x4*cos(x3)/(1-x2*parameters(1));
%             x2d=x4*sin(x3);
%             x3d=x4/parameters(2)*tan(x5/16)-parameters(1)*x4*cos(x3)/(1-x2*parameters(1));
%             x4d=parameters(3)*(u1-x4);
%             x5d=parameters(4)*(u2-x5);
%             
%             x1N=parameters(5)*t;
%             x2N=0;
%             x3N=0;
%             x4N=parameters(5);
%             x5N=16*atan(parameters(1)*parameters(2)); % we want constant curvature
%             
%             u1N=parameters(5);
%             u2N=16*atan(parameters(1)*parameters(2));
%             
%             A_syms=jacobian([x1d,x2d,x3d,x4d,x5d],[x1 x2 x3 x4 x5]);
%             B_syms=jacobian([x1d,x2d,x3d,x4d,x5d],[u1 u2]);
%             
%             A=double(subs(A_syms,[x1 x2 x3 x4 x5],[x1N x2N x3N x4N x5N]));
%             B=double(subs(B_syms,[u1 u2],[u1N u2N]));
%-----------------------------------------------------------------------

            A = double([0,parameters(1)*parameters(5),0,1,0;...
                 0,0,parameters(5),0,0;...
                 0,(-parameters(1)^2)*parameters(5),0,0,(parameters(5)*(((parameters(2)^2)*parameters(1)^2)/16 + 1/16))/parameters(2);...
                 0,0,0,-parameters(3),0;...
                 0,0,0,0,-parameters(4)]);
            B = double([0,0;0,0;0,0;parameters(3),0;0,parameters(4)]);

            C = double(eye(size(A)));
            D = double([0,0;0,0;0,0;0,0;0,0]);
            
            varargout = {A,B,C,D};
        end        
        %
        function varargout = getDiscreteLinearModel(A,B,C,D,sampling_time,method)
            % [PHI,GAM] =
            % getDiscreteLinearModel(A, B, C, D, SAMPLING_TIME, METHOD)
            % returns the PHI and GAMMA matrices characterizing the
            % discrete-time linear model of the system given the matrices
            % A,B,C,D of the continuous-time linear model and the desired 
            % SAMPLING_TIME.
            %
            % Additionally, the input METHOD is a string
            % indicating the method that should be used to calculate 
            % the matrices PHI and GAMMA. It can take values 
            % - Euler : Euler approximation as discretization method. 
            % - c2d : use the matlab command c2d. 
            Phi = [];
            Gamma = [];
            
            if strcmp(method,'Euler')
                % Calculate the discrete-time linear model using the 
                % Euler approximation.
                
                % Phi = ;
                % Gamma = ;                
                Phi = eye(size(A))+sampling_time*A;
                Gamma = sampling_time*B;
                
            elseif strcmp(method,'c2d')
                %%- Build continuous representation of the system with 'ss'
                
                 Mc = ss(A,B,C,D);
                
                %%- Calculate the discrete-time linear model of the system 
                % using the command 'c2d'
                
                 Md = c2d(Mc,sampling_time ); 
                
                %%- Extract from Md, the Phi and Gamma matrices. 
                
                 Phi = Md.A;
                 Gamma = Md.B;
                
                %%-set up output of the function 
            end
            varargout = {Phi, Gamma};
        end                
        %
        function varargout = getWorkingTrajectory(sampling_time, simulation_time, parameters)
            % [NOMINAL_TRAJECTORY_X, NOMINAL_TRAJECTORY_U] =
            % getWorkingTrajectory(SAMPLING_TIME, SIMULTAION_TIME,
            % PARAMETERS)  
            % outputs the NOMINAL_TRAJECTORY_X and NOMINAL_TRAJECTORY_U
            % given the SAMPLING_TIME between data points, the
            % SIMULATION_TIME up to which the trajectory has to be created,
            % and the vector PARAMETERS with the value sof tha system's
            % parameters.
            %
            % The outputs NOMINAL_TRAJECTORY_X, and NOMINAL_TRAJECTORY_U must
            % be arrays [t | \bar{x}] and [t | \bar{u}] 
            % whose first column corresponds to the timespan of
            % the data point, and following columns store the information
            % of the states and inputs at the corresponding time.
            %
            % The defined output trajectories are meant to be imported in
            % Simulink with the "From Workspace" block. If any
            % additional doubt regarding how the data should be formated,
            % read the information provided in the mentioned simulink block.
            %
            % todos
            % - create time vector. 
            % - create the nominal states trajectory output
            % - create the control inputs nominal trajectory output
            
            %%- create time vector
            time_vector = [0:sampling_time:simulation_time].';
            
            %%-create nominal state trajectory. 
            l=ones(length(time_vector),1); % use l to get bar_x for time 0 to simulation time
            bar_x = [parameters(5)*time_vector, 0*l, 0*l, parameters(5)*l,16*atan(parameters(1)*parameters(2))*l];
            nominal_trajectory_x = [time_vector, bar_x ];
            
            %%-create nominal control input trajectory. 
            
            bar_u = [parameters(5)*l,16*atan(parameters(1)*parameters(2))*l];
            nominal_trajectory_u = [time_vector, bar_u];
            
            varargout = {nominal_trajectory_x, nominal_trajectory_u};
        end
        %
        function varargout = getInitialState(nominal_trajectory_x)
            %[X0, X0TILDE] = getInitialState(NOMINAL_TRAJECTORY_X)
            % returns the initial state X0 of the system and the
            % initial state X0TILDE of the linear models given the 
            % information on the exercise handout and the
            % NOMINAL_TRAJECTORY_X.
            %
            % The outputs should be column vectors. 
            %
            % Remember that by definition \tilde{x} = x - \overline{x}.
            
            
            %%- define the value of x0 for experiment 1
            x0_experiment_1 =[0;0;0;5;0.0029]; %[2;2;2;2;2]% 
            %x0_experiment_2 = [0;0;0;5;0.0029];
            %%- define the value of x0Tilde for experiment 1
            x0Tilde_experiment_1 = [x0_experiment_1.'-nominal_trajectory_x(1,2:end)].';
            %x0Tilde_experiment_2 = [x0_experiment_2.'-nominal_trajectory_x(1,2:end)].';
            %including the different values for different experiments as a
            %cell
            x0 = {x0_experiment_1};%,x0_experiment_2};
            x0Tilde = {x0Tilde_experiment_1};%,x0Tilde_experiment_2};
            
            %set outputs of the function 
            varargout = {x0,x0Tilde};
        end
        %
        function varargout = getOpenLoopInputSignal(sampling_time, simulation_time)
            %[INPUT_CONTROL_ACTIONS_OPEN_LOOP] = getOpenLoopInputSignal(SAMPLING_TIME, SIMULATION_TIME)
            % outputs an input sequence to be applied in open loop to the 
            % models. The desired SAMPLING_TIME between data points as
            % well as the SIMULTION_TIME are provided. 
            %
            % As in the case of GETWORKINGTRAJECTORY function, the outputs
            % are meant to be used in Simulink with "From Workspace"
            % importing modules. If any additional doubt regarding how the
            % data should be structured, read the information provuded by
            % the mentioned simulink block. 
            %
            % todo:
            % - Declare an appropriate timespan vector. 
            % - Create the input_control_actions_open_loop array with the
            % sequence of control inputs to be applied in open loop. 
            %
            %
            % Notice: alternatively, this function can output a cell with
            % several arrays showing different control sequences to be
            % applied. This would make the provided script to run the
            % simulink mdel as many times as different control sequences
            % are gathered within the cell. Meaning that several
            % experiments can be set at once. 
            %
            
            %%- Create a time vector.
            time_vector = [0:sampling_time:simulation_time].';
            
            %%- set the control sequence to be applied in open loop for the
            %%1st experiment. 
            %time_vector = [];
            l1=[1/length(time_vector):1/length(time_vector):1].';
            l=ones(length(time_vector),1);
            u1_open_loop = 5*l; %5
            u2_open_loop = 6.4*10^(-9)*l; %0.0029
            uOpenLoop_experiment_1 = [time_vector, u1_open_loop, u2_open_loop];
            
            %u1_open_loop2= 2*l1;
            %u2_open_loop2=0.5*l;
            %uOpenLoop_experiment_2 = [time_vector, u1_open_loop2, u2_open_loop2];
            
            %Include different values for the different experiments in a
            %cell.
            input_control_actions_open_loop = {uOpenLoop_experiment_1};%,uOpenLoop_experiment_2};
            
            %set output of the function
            varargout = {input_control_actions_open_loop};
        end
        %
    end
    
end

