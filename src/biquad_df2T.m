classdef  biquad_df2T < handle
    % biquad_df2T - Initialize arrays and calculate 2nd order biquad coefficient
    % filter (lowpass, highpass, bandpass, peak, lowshelf, highshelf).
    % mono_df2T method allows you to process a mono signal by blocks of samples using the
    % direct transform 2 transposed biquad structure.
    % --------------------------
    % Author:  Oscar Butler
    % Project: MBiquad
    % Date:    11.10.2023
    % --------------------------
    
    %% Properties
    properties (Access = public)
        numStages   % Biquad order
        coeffs      % Biquad coefficient buffer
        state       % Biquad buffer state containing previous output and incoming input samples
        f0          % Significant frequency
        Q           % Quality factor of the filter
        fs          % Sampling frequency
        gaindB      % Gain in dB
        type        % Filter type
        blocksize   % Buffer size
        outputBuffer% Buffer where filtered samples are stored
    end
    %% Methods
    methods
        function obj = biquad_df2T(filter_param, fs, blocksize) % Constructor
            obj.init(filter_param, fs, blocksize);
            obj.biquad_coeff_calculation;
        end
        
        function obj = init(obj, filter_param, fs, blocksize)
            obj.type = filter_param.type;
            obj.numStages = filter_param.numStages;
            obj.coeffs = zeros(1,2*obj.numStages);
            obj.state = zeros(1,2*obj.numStages);
            obj.f0 = filter_param.freqCut;
            obj.Q = filter_param.Q;
            obj.fs = fs;
            obj.gaindB = filter_param.gaindB;
            obj.blocksize = blocksize;
            obj.outputBuffer = zeros(1,blocksize);
        end
        
        function obj = biquad_coeff_calculation(obj)
            w0 = 2*pi*obj.f0/obj.fs;
            alpha = sin(w0)/(2*obj.Q);
            switch (obj.type)
                case 0 % LP
                    A = 10^(obj.gaindB/20.0);
                    b0 = A*(1.0 - cos(w0))/2.0;
                    b1 = A*(1.0 - cos(w0));
                    b2 = A*(1.0 - cos(w0))/2.0;
                    a0 = 1.0 + alpha;
                    a1 = -2.0*cos(w0);
                    a2 = 1.0 - alpha;
                case 1 % HP
                    A = 10^(obj.gaindB/20.0);
                    b0 = A*(1.0 + cos(w0))/2.0;
                    b1 = -A*(1.0 + cos(w0));
                    b2 = A*(1.0 + cos(w0))/2.0;
                    a0 = 1.0 + alpha;
                    a1 = -2.0*cos(w0);
                    a2 = 1.0 - alpha;
                case 2 % BP
                    A = 10^(obj.gaindB/20.0);
                    b0 = A*alpha;
                    b1 = 0.0;
                    b2 = -A*alpha;
                    a0 = 1.0 + alpha;
                    a1 = -2.0*cos(w0);
                    a2 = 1.0 - alpha;
                case 3 % PK
                    A = 10^(obj.gaindB/40.0);
                    b0 = 1.0 + alpha*A;
                    b1 = -2.0*cos(w0);
                    b2 = 1.0 - alpha*A;
                    a0 = 1.0 + alpha/A;
                    a1 = -2.0*cos(w0);
                    a2 = 1.0 - alpha/A;
                case 4 % LS
                    A = 10^(obj.gaindB/40.0);
                    b0 = A*((A+1.0) - (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha);
                    b1 = 2.0*A*((A-1.0) - (A+1.0)*cos(w0));
                    b2 = A*((A+1.0) - (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha);
                    a0 = (A+1.0) + (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha;
                    a1 = -2.0*((A-1.0) + (A+1.0)*cos(w0));
                    a2 = (A+1.0) + (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha;
                case 5 % HS
                    A = 10^(obj.gaindB/40.0);
                    b0 = A*((A+1.0) + (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha);
                    b1 = -2.0*A*((A-1.0) + (A+1.0)*cos(w0));
                    b2 = A*( (A+1.0) + (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha);
                    a0 = (A+1.0) - (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha;
                    a1 = 2.0*((A-1.0) - (A+1.0)*cos(w0));
                    a2 = (A+1.0) - (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha;
                otherwise
                    b0 = 1.0;
                    b1 = 0.0;
                    b2 = 0.0;
                    a0 = 1.0;
                    a1 = 0.0;
                    a2 = 0.0;
            end
            b0 = b0/a0;
            b1 = b1/a0;
            b2 = b2/a0;
            a1 = a1/a0;
            a2 = a2/a0;
            
            obj.coeffs(1) = single(b0);
            obj.coeffs(2) = single(b1);
            obj.coeffs(3) = single(b2);
            obj.coeffs(4) = single(a1);
            obj.coeffs(5) = single(a2);
        end
        
        function obj = mono_df2T(obj, inputBuffer)
            % Load coefficients
            b0 = obj.coeffs(1);
            b1 = obj.coeffs(2);
            b2 = obj.coeffs(3);
            a1 = obj.coeffs(4);
            a2 = obj.coeffs(5);
            
            for stage=1:2:obj.numStages*2
                d1 = obj.state(stage);
                d2 = obj.state(stage+1);
                for i=1:obj.blocksize
                    inSample = inputBuffer(i);
                    
                    outSample = b0 * inSample + d1;
                    d1 = b1 * inSample - a1 * outSample + d2;
                    d2 = b2 * inSample - a2 * outSample;
                    
                    obj.outputBuffer(i) = outSample;
                end
                
                % Save new state variable by incrementing state pointer
                obj.state(stage) = d1;
                obj.state(stage+1) = d2;
                
                % The current stage input is given as the output to the next stage
                inputBuffer = obj.outputBuffer;
            end
        end
    end
end