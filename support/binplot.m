function [l1,l2]=binplot(V,x,col,type,dim)

if strcmp(type,'Area')
    bins=10.^[6:0.25:10];
    midpoints=10.^[6.125:0.25:9.875];
    scale='log';
elseif strcmp(type,'Elev')
    bins=[200 600:200:2400 3000];
    midpoints=[400 700:200:2300 2700];
    scale='linear';
elseif strcmp(type,'Dist')
    bins=[0:1e4:1.5e5];
    midpoints=[5e3:1e4:1.45e5];
    scale='linear';
end
    
    Xdistr=zeros(length(midpoints),dim);
    for i=1:length(midpoints)
        tmp=find(V>=bins(i) & V<bins(i+1));
        xtmp=x(tmp);
        if dim==3
            Xdistr(i,:)=quantile(xtmp,[0.025 0.5 0.975]);
        elseif dim==5
        Xdistr(i,:)=quantile(xtmp,[0.025 0.25 0.5 0.75 0.975]);
        end
    end
    Xdistr(Xdistr==0)=1e-12;
    
    if dim==3
        l1=patch([midpoints flip(midpoints)],[Xdistr(:,1); flip(Xdistr(:,3))],col,'FaceAlpha',0.25,'edgecolor','n');
    set(gca,'xscale',scale,'tickdir','out'); hold on
    l2=plot(midpoints,Xdistr(:,2),'color',col);
    elseif dim==5
    l1=patch([midpoints flip(midpoints)],[Xdistr(:,1); flip(Xdistr(:,5))],col,'FaceAlpha',0.25,'edgecolor','n');
    set(gca,'xscale',scale,'tickdir','out'); hold on
    patch([midpoints flip(midpoints)],[Xdistr(:,2); flip(Xdistr(:,4))],col,'FaceAlpha',0.4,'edgecolor','n');
    l2=plot(midpoints,Xdistr(:,3),'color',col);
    end
end