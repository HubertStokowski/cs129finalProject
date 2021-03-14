function F = fitlorentzian(x,y,p0,varargin)
% F = fitlorentzian(xdata,ydata,p0,[polyterms])
% fits a lorentzian to a normalized transmission spectrum
% returns F.pp = [amplitude center width]
% and F.yfit, the lorentzian over the domain with parameters P
% parameters0 = [amplitude0 center0 width0]
% leave p0 blank to autofit

c0 = 299792458; %m/s

[x ind] = sort(x(:));
y = y(ind);
y = y(:);

% if p0 = [], guess parameters
if (length(p0)==0)
    num = 50; % set to 50 in original code and wasn't working.		    meanlength = ceil(length(y)/50);
    meanlength = ceil(length(y)/num);
    ysmooth = smooth(y,meanlength);
    ddysmooth = smooth(diff(smooth(diff(ysmooth),meanlength)),meanlength);
    ddysmooth = ddysmooth(meanlength:end-meanlength);
    if(abs(min(ddysmooth))<abs(max(ddysmooth))) % lorentzian is pointing down
        p0(1) = -abs(min(ysmooth)-max(ysmooth));
        p0(2) = x(find(ysmooth==min(ysmooth),1));
    else % lorentzian is pointing up
        p0(1) = abs(min(ysmooth)-max(ysmooth));
        p0(2) = x(find(ysmooth==max(ysmooth),1));
    end
    segright = ysmooth(find(x==p0(2)):end);
    segright = segright - mean([max(segright) min(segright)]);
    segleft = ysmooth(1:find(x==p0(2)));
    segleft = segleft - mean([max(segleft) min(segleft)]);
    xcrossright = x(find(sign(segright)~=sign(segright(1)))+length(segleft)-1);
    xcrossleft = x(find(sign(segleft)~=sign(segleft(end))));
    p0(3) = diff([xcrossleft(end) xcrossright(1)]);
end

% center and rescale data so x and y range from -1 - 1 always
Tx = [ (max(x)+min(x))/2 (max(x)-min(x))/2  ];
Ty = [ (max(y)+min(y))/2 (max(y)-min(y))/2 ];
x = (x-Tx(1))/Tx(2);
y = (y-Ty(1))/Ty(2);
p0(1) = p0(1)/Ty(2);
p0(2) = (p0(2)-Tx(1))/Tx(2);
p0(3) = p0(3)/Tx(2);

opt = optimset('maxfunevals',10000,'maxiter',10000,'display','off');
pp = lsqcurvefit(@(pp,xdata) polylorentzian(pp,xdata,y,varargin{:}),p0,x,y,[],[],opt);

% uncenter and rescale data
[yfit bgfit] = polylorentzian(pp,x,y,varargin{:});
x = x*Tx(2)+Tx(1);
y = y*Ty(2)+Ty(1);
pp(1) = pp(1)*Ty(2);
pp(2) = pp(2)*Tx(2)+Tx(1);
pp(3) = abs(pp(3)*Tx(2));
yfit = yfit*Ty(2)+Ty(1);
bgfit = bgfit*Ty(2)+Ty(1);

% parameters
F.pp = pp;
F.yfit = yfit;
F.bgfit = bgfit;
F.x = x;
F.y = y;

function [ y bg ] = polylorentzian(p0,xdata,ydata,varargin)
if(~isempty(varargin)&&varargin{1}>=0)
    if(varargin{1}<=1)
        pp = zeros(1,varargin{1}+1);
        pp(end) = (ydata(1)+ydata(end))/2;
        opt = optimset('maxfunevals',10000,'maxiter',10000,'display','off');
        pp = lsqcurvefit(@(pp,xdata) polyval(pp,xdata) + p0(1)./(1+((xdata-p0(2))/(p0(3)/2)).^2),pp,xdata,ydata,[],[],opt);
        bg = polyval(pp,xdata);
    else
        stdevs = 5;
        Nbg = [1:find(xdata>p0(2)-stdevs*p0(3),1) find(xdata>p0(2)+stdevs*p0(3),1):length(xdata)];
        pp = polyfit(xdata(Nbg),ydata(Nbg),varargin{1});
        bg = polyval(pp,xdata);
    end
    y = bg + p0(1)./(1+((xdata-p0(2))/(p0(3)/2)).^2);
else
    bg = zeros(size(ydata));
    y = p0(1)./(1+((xdata-p0(2))/(p0(3)/2)).^2);
end