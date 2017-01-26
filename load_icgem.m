function [out,header,trnd,acos,asin]=load_icgem(filename)
%This function is an adaptation of icgem2mat.m from rotating_3d_globe, by
%Ales Bezdek, which can be found at:
%
%http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/
%
%The original header is transcribed below.
%
%J.Encarnacao (j.g.deteixeiradaencarnacao@tudelft.nl) 11/2013
%
% ICGEM2MAT   Reads geopotential coefficients from an ICGEM file and saves them in a mat file.
%
% Usage:
%
%       icgem2mat
%
% finds all the ICGEM files (*.gfc) in the current directory,
% reads the geopotential coefficients, transforms them into Matlab variables:
%       header...structure with Icgem header information
%       out.C(n+1,m+1), out.S(n+1,m+1)...harmonic coefficients C(n,m), S(n,m)
%
% The new mat file with the same name is moved into 'data_icgem' subdirectory;
% the original gfc file is moved into 'data_icgem/gfc/' subdirectory.
%
% Add the 'data_icgem' folder into your Matlab path.
% The model coefficients are then loaded by typing, e.g.:
%
%          load egm2008
%
% To display the C(2,0) zonal term type
%
%          out.C(3,1)
%
%
% See also compute_geopot_grids
%
% Ales Bezdek, bezdek@asu.cas.cz, 11/2012
%
% clear
% NMAX=360;
% NMAX=1e100;  %it is possible to limit the maximum degree read from the gfc file
% adr_data='./';
% adr_kam='./data_icgem/';
%
% seznam_soub=dir(adr_data);
% soub={seznam_soub.name};   %cell with filenames
% for i=1:length(soub)
%    jm=soub{i};
%    if length(jm)>4 && strcmpi(jm(end-3:end),'.gfc')
%       soub1=jm(1:end-4);
%       fprintf('Gfc file processed: %s\n',soub1);
%       filename=[adr_data soub1 '.gfc'];

  %open file
  fid=fopen(filename);
  % init header
  header=struct(...
      'product_type',           '',...
      'modelname',              '',...
      'model_content',          '',...
      'earth_gravity_constant', [],...
      'radius',                 [],...
      'max_degree',             [],...
      'errors',                 '',...
      'norm',                   '',...
      'tide_system',            '',...
      'filename',               filename...
  );
  % Read header
  s=fgets(fid);
  fn=fieldnames(header);
  while(strncmp(s, 'end_of_head', 11) == 0 && sum(s)>=0)
    for j=1:numel(fn)
      f=fn{j};
      if (keyword_search(s,f))
        valuestr=strtrim(s(length(f)+1:end));
        switch class(header.(f))
          case 'double'
            header.(f)=str2double(strrep(valuestr,'D','e'));
          case 'char'
            header.(f)=valuestr;
          otherwise
          error([mfilename,': cannot handle class ',class(header.(f)),'.'])
        end
      end
    end
    s=fgets(fid);
  end
  if s<0
    error([mfilename,'Problem reading the gfc file.'])
  end
  % sanity on max degree
  if isempty(header.max_degree)
    error([mfilename,': could not determine maximum degree of model ''',filename,'''.'])
  end
  % make room for coefficients
  out=struct('C',zeros(header.max_degree+1),'S',zeros(header.max_degree+1),'t0',[]);
  ei=struct('C',zeros(header.max_degree+1),'S',zeros(header.max_degree+1));
  trnd=struct('C',[],'S',[]);
  acos=struct('C',[],'S',[]);
  asin=struct('C',[],'S',[]);
  %iterators
  i_t0=0;
  i_trnd=0; %pocet clenu s trendem
  i_acos=0; %pocet clenu
  i_asin=0; %pocet clenu
  i_gfc=0;
  %read data
  s=fgets(fid);
  while (s>=0)
    %skip empty lines
    if numel(s)<5
      s=fgets(fid);
      continue
    end
    %retrieve keywords, degree and order
    x=str2num(s(5:end));
    n=x(1)+1;
    m=x(2)+1;
    if strcmp(s(1:4),'gfct')
      i_t0=i_t0+1;
      out.C(n,m)=x(3);
      out.S(n,m)=x(4);
      [yr,mn,dy]=ymd2cal(x(end)/1e4);
      yrd=jd2yr(cal2jd(yr,mn,dy));
      if isempty(out.t0)
        out.t0=zeros(grep_nr_occurences(filename,'gfct'));
        if numel(out.t0)==0; error([mfilename,'Problem with gfct']); end
      end
      out.t0(i_t0,:)=[n m yrd]; 
      if (strcmp(header.errors,'formal') || ...
          strcmp(header.errors,'calibrated') || ...
          strcmp(header.errors,'calibrated_and_formal'))
        ei.C(n,m)=x(5);
        ei.C(n,m)=x(6);
      end
    elseif strcmp(s(1:3),'gfc')
      out.C(n,m)=x(3);
      out.S(n,m)=x(4);
      if (strcmp(header.errors,'formal') || ...
          strcmp(header.errors,'calibrated') || ...
          strcmp(header.errors,'calibrated_and_formal'))
        ei.C(n,m)=x(5);
        ei.C(n,m)=x(6);
      end
      i_gfc=i_gfc+1;
    elseif strcmp(s(1:4),'trnd') || strcmp(s(1:3),'dot')
      if isempty(trnd.C)
        trnd.C=zeros(grep_nr_occurences(filename,'trnd')+grep_nr_occurences(filename,'dot'),3); trnd.S=trnd.C;
        if numel(trnd.C)==0; error([mfilename,'Problem with trnd']); end
      end
      i_trnd=i_trnd+1;
      trnd.C(i_trnd,:)=[n m x(3)]; 
      trnd.S(i_trnd,:)=[n m x(4)]; 
    elseif strcmp(s(1:4),'acos')
      if isempty(acos.C)
        acos.C=zeros(grep_nr_occurences(filename,'acos'),4); acos.S=acos.C;
        if numel(asin.C)==0; error([mfilename,'Problem with acos']); end
      end
      i_acos=i_acos+1;
      acos.C(i_acos,:)=[n m x(3) x(end)];
      acos.S(i_acos,:)=[n m x(4) x(end)];
    elseif strcmp(s(1:4),'asin')
      if isempty(asin.C)
        asin.C=zeros(grep_nr_occurences(filename,'asin'),4); asin.S=asin.C;
        if numel(asin.C)==0; error([mfilename,'Problem with asin']); end
      end
      i_asin=i_asin+1;
      asin.C(i_asin,:)=[n m x(3) x(end)];
      asin.S(i_asin,:)=[n m x(4) x(end)];
    else
      error([mfilename,'A problem occured in gfc data.']);
    end
    s=fgets(fid);
  end
  fclose(fid);
  %handle the tide system
  switch header.tide_system
    case 'zero_tide'
      %do nothing, this is the default
    case {'free_tide','tide_free'}
      out.C(3,1)=out.C(3,1)-4.173e-9;
      header.tide_system='zero_tide';
    case 'mean_tide'
      %http://mitgcm.org/~mlosch/geoidcookbook/node9.html
      out.C(3,1)=out.C(3,1)+1.39e-8;
      header.tide_system='zero_tide';
    otherwise
      %the tide system is not documented, so make some assumptions
      switch header.modelname
        case 'GROOPS'
          %do nothing, Norber Zehentner's solutions are zero tide
        otherwise
          error([mfilename,': unknown tide system ''',header.tide_system,'''.'])
      end
  end
end

function out=keyword_search(line,keyword)
    out=strncmp(strtrim(line),       keyword,         length(keyword)) || ...
        strncmp(strtrim(line),strrep(keyword,' ','_'),length(keyword));
end
