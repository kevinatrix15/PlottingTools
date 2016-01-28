clear all;
clc;
clf;
csvfiles = dir('*_Temperature.csv');
[~,index] = sortrows([csvfiles.datenum].'); csvfiles = csvfiles(index);
numfiles = length(csvfiles);
% numfiles = 10;

for k=1:numfiles
    Data=csvread(csvfiles(k).name);
    p=Data(:,1:4);
    xyrange=find(p(:,3)==min(p(:,3)));
    xylength=size(xyrange,1);
    zrange=unique(p(:,3));
    zlength=size(zrange,1);
    count=0;
    liquidus_T=1923;
    mean_T=[];
    pos=[];
    f=[];
    %===============================================================
    % average temperature for each z where t is less than liquidus_T
    %===============================================================
    [e,f]=find(p((zlength-1)*xylength+1:zlength*xylength,4)>1923);
    [c,d]=find(p((zlength-1)*xylength+1:zlength*xylength,4)-1923>0 & p((zlength-1)*xylength+1:zlength*xylength,4)-1923<=50);
    if ~isempty(c)==1
        subplot(2,2,1)
        plot(p(e,1),p(e,2),'*');
        hold on
        plot(p(c,1),p(c,2),'o');
        hold off
        subplot(2,2,2)
        plot3(p((zlength-1)*xylength+1:zlength*xylength,1),p((zlength-1)*xylength+1:zlength*xylength,2),p((zlength-1)*xylength+1:zlength*xylength,4))
        for i=1:size(c,1)
            xmark=p(c(1),1);
            ymark=p(c(1),2);
            [a,b]=find(p(:,1)==xmark & p(:,2)==ymark);
            pos(:,i)=a;
            count=count+1;
        end
        for i=1:zlength
            mean_T(i,1)=mean(p(pos(i,:),4));
        end
        % add more points in the coarse region in order to ease curve fit
        mintol=1E-4;
        tol=2E-5;
        zdiff=diff(zrange);
        Tdiff=diff(mean_T);
        [g,h]=find(zdiff>=mintol);
        if ~isempty(g)
            index=0;
            for tt=1:length(g)
                zintpts(tt)=floor(zdiff(g(tt))/tol);
                zintpos=[];
                Tintpos=[];
                for pp=1:zintpts(tt)
                    zintpos(pp,1)=zrange(index+g(tt))+pp*tol;
                    Tintpos(pp,1)=mean_T(index+g(tt))+pp*tol*Tdiff(g(tt))/zdiff(g(tt));
                end
                zrange=[zrange(1:g(tt)+index);zintpos;zrange(g(tt)+index+1:end)];
                mean_T=[mean_T(1:g(tt)+index);Tintpos;mean_T(g(tt)+index+1:end)];
                index=index+zintpts(tt);
            end
        end
        % zrange=sort(zrange);
        zrange=-zrange+max(zrange);
        mean_T=(mean_T-min(mean_T))/(max(mean_T)-min(mean_T));
        % cftool;
        f=fit(zrange,mean_T,'exp(-b*x)');
        bcoef(k,1)=coeffvalues(f);
        confrange(k,:)=confint(f);
        

        figure(4)
        plot(mean_T,-zrange,'ro');
        hold on;
        bcoef(~bcoef)=[];
        syms x;
        curve_T = eval(subs('exp(-bcoef*x)', x, zrange));
        plot(curve_T,-zrange, '--k', 'LineWidth', 2.25);
        xlabel('$$\frac{T(z) - T_h}{T_0 - T_h}$$','Interpreter','latex', 'FontSize', 14);
        ylabel('$$z [m]$$','Interpreter','latex', 'FontSize', 14);
        legend('T_{sim}', 'h(z)', 'Interpreter', 'latex')
        axis([-0.05 1 -1e-3 0])
        hold off;
        pause;
        
    end
end

mean_b=mean(bcoef);


