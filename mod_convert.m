function out=mod_convert(in,from,to,GM,R,factor)
% OUT=MOD_CONVERT(IN,FROM,TO) converts model IN from scale FROM into scale
% TO.
%
%   MOD_CONVERT(IN,FROM,TO,FACTOR) additionally, the output is scalled by
%   FACTOR. This is useful to request scales in units different than SI,
%   i.e. milimeters geoid height instead of meters.
%
%   Supported units are the following,
%       'non-dim'   - non-dimensional Stoked coefficients.
%       'eqwh'      - equivalent water height [m]
%       'geoid'     - geoid height [m]
%       'potential' - [m2/s2]
%       'gravity'   - [m /s2], if input represents the disturbing potential,
%                     then the output represents the gravity disturbances.
%                     Otherwise, it represents the gravity accelerations
%                     (-dU/dr).
%       'anomalies' - If input represents the disturbing potential, then
%                     the output represents the gravity anomalies.
%                     Otherwise, it represents (-dU/dr - 2/r*U).
%       'vert-grav-grad' - vertical gravity gradient.

if ~exist('GM','var') || isempty(GM)
    GM=398600441500000; % m^3 s^-2
end
if ~exist('R','var') || isempty(R)
    R =6378136.460;     % Earth's equatorial radius [m]
end
if ~exist('factor','var') || isempty(factor)
    factor=1;
end

%getting scale
scale=mod_convert_aux(from,to,size(in.C,1),GM,R)*factor;

%scaling
out=struct(...
  'C',in.C.*scale,...
  'S',in.S.*scale...
);

function scale=mod_convert_aux(from,to,N,GM,R)

%universal grav const
G=6.6732e-11;       %m3/kg/s2

R_av=6371000;       % Earth's mean radius [m]
rho_earth= (GM/G) / (4/3*pi*R_av^3); %kg/m3 (average density of the Earth)
rho_water=1000;

%defining love numbers
love=[  0       0.000;...
        1       0.027;...
        2      -0.303;...
        3      -0.194;...
        4      -0.132;...
        5      -0.104;...
        6      -0.089;...
        7      -0.081;...
        8      -0.076;...
        9      -0.072;...
        10     -0.069;...
        12     -0.064;...
        15     -0.058;...
        20     -0.051;...
        30     -0.040;...
        40     -0.033;...
        50     -0.027;...
        70     -0.020;...
        100    -0.014;...
        150    -0.010;...
        200    -0.007];

switch lower(from)
    case 'non-dim'
        %no scaling
        pre_scale=ones(N);
    otherwise
        %need to bring these coefficients down to 'non-dim' scale
        pre_scale = mod_convert_aux('non-dim',from,N,GM,R);
end

switch lower(to)
    case 'non-dim'
        %no scaling
        pos_scale=ones(N);
    case 'eqwh' %[m]
        pos_scale=zeros(N);
        %converting Stokes coefficients from non-dimensional to equivalent water layer thickness
        for i=1:N
            deg=i-1;
            lv=interp1(love(:,1),love(:,2),deg,'linear','extrap');
            pos_scale(i,:)=R*(rho_earth/rho_water) * 1/3 * (2*deg+1)/(1+lv);
        end
    case 'geoid' %[m]
        pos_scale=ones(N)*R;
    case 'potential' %[m2/s2]
        pos_scale=ones(N)*GM/R;
    case 'gravity' %[m/s2]
        %If input represents the disturbing potential, then the output
        %represents the gravity disturbances. Otherwise, it represents the
        %gravity accelerations (-dU/dr).
        deg=(0:N-1)'*ones(1,N);
        pos_scale=-GM/R^2*(deg+1);
    case 'anomalies'
        %If input represents the disturbing potential, then the output
        %represents the gravity anomalies. Otherwise, it represents:
        %-dU/dr - 2/r*U
        deg=(0:N-1)'*ones(1,N);
        pos_scale=GM/R^2*max(deg-1,ones(size(deg)));
    case 'vert-grav-grad'
        deg=(0:N-1)'*ones(1,N);
        pos_scale=GM/R^3*(deg+1).*(deg+2);
    otherwise
        error([mfilename,': unknown scale ',to])
end

scale=pos_scale./pre_scale;

