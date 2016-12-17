function   plot_pendulum(t,y,l)
figure(2)
clf
hold
plot(l*sin(y),-l*cos(y),'ko','Markersize',12,'MarkerFaceColor',[0 0 0])
line([0 l*sin(y)],[0 -l*cos(y)],'Linewidth', 2,'Color',[0 0 0])
axis([-12 12 -24 0])

end

