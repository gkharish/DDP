clear all
close all

X = [0.0;0.0;0.0*1e5;2.0*1e5]
X2 = [0.0;0.0;1.0;4.0*1e5]
U1 = 2.0
U2 = 2.0
U1dot = 0.0
U2dot = 0.0
Xlist = [];
X2list = [];
final_time = 10;
dt = 5e-3;
N = final_time/dt;


for i=1:N
  Xdot = xdotpneu(X,U1,U2);
  % Xddot = computeXddot(X,Xdot,U,Udot);
  X = X + dt*Xdot;
  % X2 = X2 + dt*Xdot + ((dt^2)/2)*Xddot;
  Xlist = [Xlist, X];
  % X2list = [X2list, X2];
end


%xode = lsode ("xdotpneu", X, (t = linspace(0, final_time, N)));

figure()
subplot(221)
plot(linspace(0,final_time,N),Xlist(1,:))%,linspace(0,final_time,N),X2list(1,:))
subplot(222)
plot(linspace(0,final_time,N),Xlist(2,:))%,linspace(0,final_time,N),X2list(2,:))
subplot(223)
plot(linspace(0,final_time,N),Xlist(3,:))%,linspace(0,final_time,N),X2list(3,:))
subplot(224)
plot(linspace(0,final_time,N),Xlist(4,:))%,linspace(0,final_time,N),X2list(4,:))

figure()
subplot(221)
plot(t,xode(:,1))
subplot(222)
plot(t,xode(:,2))
subplot(223)
plot(t,xode(:,3))
subplot(224)
plot(t,xode(:,4))
