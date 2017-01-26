function [long,lat,grid_out]=mod_sh_synth(model,lat,NLon)
% [LONG,LAT,GRID]=MOD_SH_SYNTH(MODEL) is the low-level spherical harmonic
% synthesis routine. It should be a self-contained script. Inputs MODEL.C and
% MODEL.S are the cosine and sine coefficients, organized in lower triangle
% matrices, with constant degree in each row and constant order in each
% column, sectorial coefficients are in the diagonals, zonal coefficients
% are in the first column of matrix C (first column of S matrix is zeros).
%
%   MOD_SH_SYNTH(MODEL,LAT,NLON) with the optional inputs LAT, NLON causes
%   the synthesis to be done on the latitudes specified with vector LAT and
%   along the #NLON equidistant number of points along the longitudes
%   domain.
%
%   Output LONG is a scalar.
%   Output LAT is a vertical vector [-pi/2:pi/2].
%   Output GRID is a matrix with size [length(lat) NLon]
%
%   All angular quantities are in radians.

% Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>
% List of changes:
%
%   P.Inacio <p.m.g.inacio@tudelft.nl>, 10/2011, Added check for repeated
%       meridians at 0 and 2*pi. If they both exist and have the same Also
%       added a more robust check that grids are regular taking into
%       account a numerical error threshold.
%   P.Inacio <p.m.g.inacio@tudelft.nl>, 03/2012, The algorightm upon which
%       this function relies only uses the number of points along
%       longitude. Therefore, to avoid misinterpretations the input LONG
%       has been replaces with the number of points along the longitude
%       domain. Further interpolation to arbitrary sets of longitudes,
%       regular or irregular is implemented outside this routine, in
%       mod_synthesis.

  % Defaults
  %check dimensions
  if any(size(model.C) ~= size(model.S))
       error([mfilename,': inputs <model.C> and <model.S> must have the same size.'])
  end
  if size(model.C,1) ~= size(model.C,2)
      error([mfilename,': inputs <model.C> and <model.S> must be square matrices'])
  end

  %calculating max resolution
  L = size(model.C,1)-1;
  N = max([1 L]);
  % %%% OLD CODE
  % N = max([1,min([L,360])]);
  % %%%

  %building latitude domain (if not given)
  if ~exist('lat','var') || isempty(lat)
      lat=linspace(-pi/2,pi/2,N)';
  end
  %setting default number of points along latitude
  if ~exist('NLon','var') || isempty(NLon)
      NLon = 2*N;
  elseif ~isscalar(NLon) || ~isnumeric(NLon)
      error([mfilename,': input <Nlon> must be a numeric scalar.'])
  end

  %create internal longitude domain - only used for output.
  long=linspace(0,2*pi,NLon+1);
  % removing duplicate zero longitude
  long(end)=[];

  %checking (possible) inputs
  % NOTE: Requiring latitude to be vertical is a bit stupid, but it general
  %       it allows it to be distinguished from the longitude avoiding
  %       possible usage error.
  if size(lat,2) ~= 1
      error([mfilename,': input <lat> must be a vertical vector.'])
  end

  %check latitude domain
  if max(lat) > pi/2 || min(lat) < -pi/2
      error([mfilename,': input <lat> does not seem to be in radians or outside legal domain [-pi/2,pi/2].'])
  end

  %need co-latitude
  clat=pi/2-lat;

  % calculating Legendre polynomials, P{latitude}
  P=legendre_latitude(L,clat);

  %building the grids
  grid_out = zeros(length(clat),NLon);
  for i=1:length(clat)
      %calculating Fourier coefficients
      a=zeros(L+1,1);
      b=zeros(L+1,1);
      for m=0:L
          n=m:L;
          a(m+1) = sum(model.C(n+1,m+1).*P{i}(n+1,m+1));
          b(m+1) = sum(model.S(n+1,m+1).*P{i}(n+1,m+1));
      end

      %FFT
      grid_out(i,:) = real(fft(a+1i*b,NLon))';
  end
end
function out = legendre_latitude(L,lat)
  %getting legendre coefficient, per degree
  P=legendre_degree(L,lat);
  %building legendre coefficients, per latitude
  out = cell(length(lat),1);
  for i=1:length(lat)
      out{i} = zeros(L+1,L+1);
      for n=0:L
          out{i}(n+1,1:n+1) = P{n+1}(:,i)';
      end
  end
end
function out = legendre_degree(L,lat)
  if ~isvector(lat)
      error([mfilename,': input <lat> must be a vector.'])
  end
  %getting legendre coefficients, per degree
  out=cell(1,L+1);
  for n=0:L
      out{n+1}=legendre(n,cos(lat'),'norm')*2;
      out{n+1}(1,:)=out{n+1}(1,:)/sqrt(2);
  %     % P.Inacio - testing - manual normalization
  %     m = (0:n)';
  %     out{n+1}=legendre(n,cos(lat'));
  %     out{n+1}=repmat(sqrt((2*n+1)*factorial(n-m)./factorial(n+m)),[1 length(lat)]).*out{n+1};
  %     % P.Inacio - testing - no normalization
  %     out{n+1}=legendre(n,cos(lat'),'norm')*2;
  end
end