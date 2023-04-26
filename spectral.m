function [IrData,IomData]=spectral(Filename,Options);
% SPECTRAL Spectral analysis
%    [Ir,Iom] = SPECTRAL('FILENAME',OPTIONS)
%    Spectral.m - Performes a spectral analysis as developed for ecology 
%    and geosciences by Renshaw and coworkers. FILENAME is the name of a
%    grayscale JPG or TIF file (TIF is better than JPG).
%
%    Length gives the dimensions of the image in the horizontal
%
%    Options contains character flags that do the following
%    'n' No graphs are displayed
%    'b' The bin distribution is shown
%    'c' Force a recalculation
%    't' Progress is given in text at the command line
%    'p' A progress bar shows progress
%
%    For references sbout spectral analysis of 2D data see 
%    Renshaw and Ford 1984, Vegetatio 56:75-85
%    Muggulestone and Renshaw 1998, Computers & Geosciences 24
%    Couteron and Lejeune 2001, Journal of Ecology 89:616-628

% preparing the fotos to square size
% im=imread('Fig2bNew.bmp'); mnd=size(im);
% imwrite(im(1:min(mnd(1:2)),1:min(mnd(1:2)),:),'myimage.tif');

global Ir Iom Chi_up_Ang Chi_dw_Ang Chi_up_Rad Chi_dw_Rad
% matlabpool open
on=1;off=0;
Graph=on;PlotBins=off;Calc=off;Pbar=off;ProgText=off;

if nargin>1,
    Graph=1-sum(Options=='n')>0;
    PlotBins=sum(Options=='b')>0;
    Calc=sum(Options=='c')>0;    
    Pbar=sum(Options=='p')>0;    
    ProgText=sum(Options=='t')>0;        
end;

% Get Screen dimensions and set Main Window Dimensions
x = get(0,'ScreenSize');
ScreenDim=x(3:4);
MainWindowDim=floor(ScreenDim.*[0.9 0.8]);

i=length(Filename);

Ext=lower(Filename([i-2:i]));

switch Ext
    case {'jpg','tif'},
       Image=imread(Filename,Ext);            % Leest 3D variabele AM in vanuit de file, als een jpg file
       Filename=Filename([1:i-4]);
    otherwise
       Image=imread(Filename,'tif');          % Leest 3D variabele AM in vanuit de file, neemt aan een TIF file
end;

Datafile=[Filename '.mat'];

ShowSection=32;             % Size of the submatrix of I that is shown (a centre cutout)

ImInfo=imfinfo(Filename,Ext);
Colorpicture=1-strcmp(ImInfo.ColorType,'grayscale');
if Colorpicture,
   ImageBW=rgb2gray(Image);
   [m,n,d]=size(Image); % Obtaining the size of the image array;
else
   ImageBW=Image; 
   [m,n]=size(Image); 
   Image=zeros(m,n,3);
   Image(:,:,1)=ImageBW;
   Image(:,:,2)=ImageBW;
   Image(:,:,3)=ImageBW;
end;

% Obtaining the size of the image array;
[m,n]=size(ImageBW);
if m~=n,
    disp('The image file needs to be square');
    return;
end;

AngleBinSize=10;            % The size of a angle bin in the angular spectrum
Resolution= 50/n;        % 1 pixel is xxx meters;

pmax=m/2;
qmax=m/2;

Y=zeros(m,n);
X=zeros(m,n);

a=zeros(pmax,qmax*2);
b=zeros(pmax,qmax*2);
I=zeros(pmax,qmax*2);
V=zeros(pmax,qmax*2);

Y(:)=ImageBW(:);   
% The image arraw is copied to a real-array that allows calculations

X = Y - mean(Y(:)); 
% Y is rescaled, so that the average equals zero;

% Below, the periodogram is calculated

if ProgText, 
    disp('Processing image ...'); disp('Currently at   0%');
end;
if Pbar, h=waitbar(0,'Calculating periodogram');end;

[t,s]=meshgrid(1:n,1:m);    
for p=0:pmax,
   for q=0:qmax*2-1,
       a(p+1,q+1)=sum( sum(X.*cos(2*pi*(p*s/m+(q-qmax)*t/n) ) ) )/(m*n);
       b(p+1,q+1)=sum( sum(X.*sin(2*pi*(p*s/m+(q-qmax)*t/n) ) ) )/(m*n);
       I(p+1,q+1)=m*n*(a(p+1,q+1)^2+b(p+1,q+1)^2);
   end;
   if ProgText,
       disp(sprintf('\b\b\b\b\b%3.0f%%',p/pmax*100));
   end;
   if Pbar, waitbar(p/pmax,h);end;       
   drawnow;
end;
if Pbar, close(h);end;       
% The variance is calculated


V=X.*X;                         % The variance array
Vtot=sum(V(:))/(m*n);           % The total variance

Dom=AngleBinSize;               % The AngleBinSize is copied into a variable with a shorter name
om_max=180/Dom;                 % The total number of bins

[pmax qmax]=size(I);            

rmax=min([pmax qmax]);

% Variable definitions
Ir  = zeros(1,rmax);            % The bins vor the radial spectrum is defined
Irc = zeros(1,rmax);            % The array with bin sizes for the radial spectrum is defined
Iom = zeros(1,om_max);          % The bins for the angular spectrum is defined
Iomc= zeros(1,om_max);          % The array with bin sizes for the angular spectrum is defined

Chi_up_Rad = zeros(1,rmax);     % The upper 95% confidence limits, Angular spectrum
Chi_dw_Rad = zeros(1,rmax);     % The lower 95% confidence limits, Angular spectrum
Chi_up_Ang = zeros(1,om_max);   % The upper 95% confidence limits, Angular spectrum
Chi_dw_Ang = zeros(1,om_max);   % The lower 95% confidence limits, Angular spectrum

Id=zeros(pmax,qmax);

% The bins are defined spatially
[q,p]=meshgrid(1:qmax,0:pmax-1);
qs=q-0.5*qmax-1;

warning off MATLAB:divideByZero
om=ceil(atan(p./(qs))/pi*180/Dom)+om_max/2;
warning on;
om(1,(qmax/2+1):qmax)=0;

rs=round(sqrt(p.^2+(qs).^2));
rs(1,(qmax/2+1):qmax)=0;

% The angular bins are filled, and the number of cells counted
for i=1:om_max,
   Iom(i) = sum(sum((om==i).*I));
   Iomc(i)= sum(sum((om==i)));
end;    
    
% The radial bins are filled, and the number of cells counted
for i=1:rmax,
   Ir(i) = sum(sum((rs==i).*I));
   Irc(i)= sum(sum((rs==i)));
end;

% The bins are rescaled
Ir(:)=Ir(:)./Vtot./Irc(:);             % Not correct
Iom(:)=Iom(:)./Vtot./Iomc(:); 

% Calculation of the Chi-squared critical values
Chi_up_Ang(:) = 1./(2*Iomc(:)).*chi2inv(0.975,Iomc(:)*2);
Chi_dw_Ang(:) = 1./(2*Iomc(:)).*chi2inv(0.025,Iomc(:)*2);
Chi_up_Rad(:) = 1./(2*Irc(:)).*chi2inv(0.975,Irc(:)*2);             
Chi_dw_Rad(:) = 1./(2*Irc(:)).*chi2inv(0.025,Irc(:)*2);            

% Doubling the original image for viewing, the lower part is rotated upwards
Itot = zeros(qmax,qmax);
Itot(pmax:2*pmax-1,:)=I(:,:);
Itot(1:pmax,:)=fliplr(flipud(I));
Itot(1:pmax,2:qmax)=Itot(1:pmax,1:qmax-1);

% Taking a subimage from the centre of the I image, with size ShowSection
s=ShowSection;
Ishow=zeros(ShowSection,ShowSection);
Ishow(:,:)=Itot((0.5*(m-s)+1):(0.5*(m+s)),(0.5*(n-s)+1):(0.5*(n+s)));

% Finally, the figure is drawn

IomData=[Iom; Chi_up_Ang; Chi_dw_Ang];
IrData=[1:rmax; Ir; Chi_up_Rad; Chi_dw_Rad];

% matlabpool close
