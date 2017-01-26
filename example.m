
%load the coefficients
[i,header]=load_icgem('ggm05g.gfc.txt');

%convert to geoid height
g=mod_convert(i,'non-dim','geoid',header.earth_gravity_constant,header.radius);

%perform synthesis with 1 degree samplint at the equator
[long,lat,grid_out]=mod_sh_synth(g,0,360);

%plot the geoid height above the Earth mean radius of the model
plot(long*180/pi,grid_out-header.radius);