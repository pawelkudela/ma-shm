function [h1,strt,endd] = propPlotEnery(time,E,E2,thrs,trs,PRT,file)

samples = size(time,2);

h1 = figure('Name',file);
plot(time,E2/max(E2),'black','LineWidth',2); hold
plot(time,E/max(E2),'black:');
xlabel('Time [s]'); 
ylabel('Normalized signal energy'); ylim([0 1.1]); xlim([time(1),time(end)]);


[maxx ~] = max(E);
strt = 2;

while strt < samples
    if E(strt) > trs*maxx
        break
    else
     strt = strt + 1;
    end
end

endd = strt;
while endd < samples
    if E(endd) < thrs/100*maxx
        break
    else
     endd = endd + 1;
    end
end

hold on
line([time(strt) time(strt)], [0 1],'LineWidth',1,'LineStyle','--','Color',[0.2 .2 1]);
line([time(endd) time(endd)], [0 E2(endd)/maxx],'LineWidth',1,'LineStyle','--','Color',[.2 .2 1] )
     
if PRT == 1
    print('-djpeg','-r200', [file, '.jpg']);
end


