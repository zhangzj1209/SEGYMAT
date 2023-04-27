function [s,t] = sweep(fmin,fmax,dt,tmax,staper,etaper,ltaper,amplitude,theta0)
% SWEEP: generate a linear Vibroseis sweep
%
% [s,t]= sweep(fmin,fmax,dt,tmax,staper,etaper,ltaper,amplitude,theta0)
% [s,t]= sweep(); %use all defaults
% [s,t]= sweep([],[],2e-3)%use defaults, but override sample rate
% NOTE that ltaper=2.5 better matches the high frequency end of a sweep 
%      with the same parameters generated by a Pelton Vibe controller
%
% SWEEP generates a linear synthetic Vibroseis sweep for the 
% specified passband. Reference: Aldrige, D.F., 1992, Mathematics of linear sweeps, vol 28(1),
% 62-68, http://csegjournal.com/assets/pdfs/archives/1992_06/1992_06_math_linear_sweeps.pdf
%
% fmin= minimum swept frequency in Hz
%       default = 10.;
% fmax= maximum swept frequency in Hz
%       default = 100.;
% dt= time sample rate in seconds
%       default = 1e-3 (1 ms)
% tmax= sweep length in seconds
%       default = 10.;
% staper= length of start taper (cos) in seconds
%       default = 0.2
% etaper= length of end taper (cos) in seconds
%       default = staper
% ltaper= percentage reduction in amplitude over length of sweep
%       default = 0.
% amplitude= maximum amplitude of the sweep
%       default = 1.
% theta0= start phase of the sweep
%       default = -pi/2.
%
% s= output linear sweep
% t= output time coordinate vector 
%
% by K.W. Hall, Oct 2020. Original function is now called sweep_old.m
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

% check inputs
%function [s,t] = ksweep(fmin,fmax,dt,tmax,staper,etaper,ltaper,amplitude,theta0)

if nargin<1 || isempty(fmin)
    fmin = 10.;
end
if nargin<2 || isempty(fmax)
    fmax = 100.;
end
if nargin<3 || isempty(dt)
    dt = 1e-3;
end
if nargin<4 || isempty(tmax)
    tmax = 10.;
end
if nargin<5 || isempty(staper)
    staper=0.2;
end
if nargin<6 || isempty(etaper)
    etaper=staper;
end
if nargin<7 || isempty(ltaper)
    ltaper=0.;
end
if nargin<8 || isempty(amplitude)
    amplitude=1.;
end
if nargin<9 || isempty(theta0)
    theta0 = -pi/2; 
end

%% create time vector
t = (0:tmax/dt)*dt;

%% create amplitude function
a = amplitude*ones(size(t));

%start taper
if staper>0
    a(t<=staper) = (amplitude/2) * (1 - cos(pi*t(t<=staper)/staper));
end
%end taper
if etaper>0
    a(t>=tmax-etaper) = (amplitude/2) * (1 + cos(pi*(t(t>=tmax-etaper)-tmax+etaper)/etaper));
end

%linear taper
a = a.*linspace(1,1-ltaper/100,length(t));

%% create sweep
s = a.*cos(theta0 +2*pi*fmin*t +pi*((fmax-fmin)/tmax)*t.^2);



