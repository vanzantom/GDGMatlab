function [x,y]=square_foursub(bs,s)
%SQUAREG   Gives geometry data for the squareg PDE model
%   The square [0,1]x[0,1] is divided into 4 subdomains.
%   NE=square_foursub gives the number of boundary segment
%
%   D=square_foursub(BS) gives a matrix with one column for each boundary segment
%   specified in BS.
%   Row 1 contains the start parameter value.
%   Row 2 contains the end parameter value.
%   Row 3 contains the number of the left hand region.
%   Row 4 contains the number of the right hand region.
%
%   [X,Y]=square_foursub(BS,S) gives coordinates of boundary points. BS specifies the
%   boundary segments and S the corresponding parameter values. BS may be
%   a scalar.

% Copyright 1994-2016 The MathWorks, Inc.

nbs=12;

if nargin==0
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 0 0 0 0 0 0 0 0   % start parameter value
  1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2  % end parameter value
  0 2 4 0 3 0 0 0 0 3 0 0 % left hand region
  1 1 1 1 4 4 4 2 2 2 3 3 % right hand region
];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs)
  error(message('pde:squareg:InvalidBs'))
end

if nargin==1
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1
  bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2)
  error(message('pde:squareg:SizeBs'));
end

if ~isempty(s)

% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=interp1([d(1,1),d(2,1)],[0 1/2],s(ii));
y(ii)=interp1([d(1,1),d(2,1)],[1 1],s(ii));
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=interp1([d(1,2),d(2,2)],[1/2 1/2],s(ii));
y(ii)=interp1([d(1,2),d(2,2)],[1 1/2],s(ii));
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=interp1([d(1,3),d(2,3)],[1/2 0],s(ii));
y(ii)=interp1([d(1,3),d(2,3)],[1/2 1/2],s(ii));
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=interp1([d(1,4),d(2,4)],[0 0],s(ii));
y(ii)=interp1([d(1,4),d(2,4)],[1/2 1],s(ii));
end
% boundary segment 5
ii=find(bs==5);
if length(ii)
x(ii)=interp1([d(1,5),d(2,5)],[1/2 1/2],s(ii));
y(ii)=interp1([d(1,5),d(2,5)],[1/2 0],s(ii));
end
% boundary segment 6
ii=find(bs==6);
if length(ii)
x(ii)=interp1([d(1,6),d(2,6)],[1/2 0],s(ii));
y(ii)=interp1([d(1,6),d(2,6)],[0 0],s(ii));
end
% boundary segment 7
ii=find(bs==7);
if length(ii)
x(ii)=interp1([d(1,7),d(2,7)],[0 0],s(ii));
y(ii)=interp1([d(1,7),d(2,7)],[0 1/2],s(ii));
end
% boundary segment 8
ii=find(bs==8);
if length(ii)
x(ii)=interp1([d(1,8),d(2,8)],[1/2 1],s(ii));
y(ii)=interp1([d(1,8),d(2,8)],[1 1],s(ii));
end
% boundary segment 9
ii=find(bs==9);
if length(ii)
x(ii)=interp1([d(1,9),d(2,9)],[1 1],s(ii));
y(ii)=interp1([d(1,9),d(2,9)],[1 1/2],s(ii));
end
% boundary segment 10
ii=find(bs==10);
if length(ii)
x(ii)=interp1([d(1,10),d(2,10)],[1 1/2],s(ii));
y(ii)=interp1([d(1,10),d(2,10)],[1/2 1/2],s(ii));
end
% boundary segment 11
ii=find(bs==11);
if length(ii)
x(ii)=interp1([d(1,11),d(2,11)],[1 1],s(ii));
y(ii)=interp1([d(1,11),d(2,11)],[1/2 0],s(ii));
end
% boundary segment 12
ii=find(bs==12);
if length(ii)
x(ii)=interp1([d(1,12),d(2,12)],[1 1/2],s(ii));
y(ii)=interp1([d(1,12),d(2,12)],[0 0],s(ii));
end
end