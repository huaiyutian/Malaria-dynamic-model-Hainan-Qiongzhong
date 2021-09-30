function ss = Malaria2(theta,data)

time   = data.ydata(:,1);
ydata  = data.ydata(:,2);
ydata1  = data.ydata(1:15,3);
xdata  = data.xdata;

y0 = theta(end-1:end);
ymodel = Malaria3(time,theta,y0,xdata);

ss =sum((sqrt(ymodel(:,1)) - sqrt(ydata(:,1))).^2);%...
