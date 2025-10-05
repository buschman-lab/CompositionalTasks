function th = AddSig(pval,pos,numcomp,nsoffset,show_p,angle)
    if nargin <3
        numcomp = 1;
    end
    if nargin<4
        nsoffset = 0.01;
    end
    if nargin <5
        show_p = 0;
    end
    if nargin <6
        angle = 0;
    end
    loc = pos(1)+(abs(pos(2)-pos(1))/2);
    if show_p        
        line([pos(1) pos(2)],[pos(3),pos(3)],'Color',[0.1 0.1 0.1],'LineWidth',2)
        if pval<=1e-16
            pval=1e-16;
            exponent = floor(log10(pval));
            str = (sprintf('p<10^{%d}',exponent));
            th = text([loc loc],[pos(3)+nsoffset pos(3)+nsoffset],str,'Color',[0.1 0.1 0.1],'FontSize',14,'FontWeight','normal','HorizontalAlignment','Center','Rotation',angle);                  
        else
            exponent = floor(log10(pval));
            if exponent<-3
                str = (sprintf('p=10^{%d}',exponent));
            else
                str = (sprintf('p=%.2g',pval));
            end
            th =text([loc loc],[pos(3)+nsoffset pos(3)+nsoffset],str,'Color',[0.1 0.1 0.1],'FontSize',14,'FontWeight','normal','HorizontalAlignment','Center','Rotation',angle);                  
        end
    else        
        line([pos(1) pos(2)],[pos(3),pos(3)],'Color',[0.1 0.1 0.1],'LineWidth',2)
        if pval<=0.001/numcomp 
            th = text([loc loc],[pos(3)+nsoffset pos(3)+nsoffset],'***','Color',[0.1 0.1 0.1],'FontSize',18,'FontWeight','normal','HorizontalAlignment','Center','Rotation',angle);
        elseif pval<=0.01/numcomp & pval>0.001/numcomp
            th = text([loc loc],[pos(3)+nsoffset pos(3)+nsoffset],'**','Color',[0.1 0.1 0.1],'FontSize',18,'FontWeight','normal','HorizontalAlignment','Center','Rotation',angle);            
        elseif pval<0.05/numcomp & pval>0.01/numcomp 
            th = text([loc loc],[pos(3)+nsoffset pos(3)+nsoffset],'*','Color',[0.1 0.1 0.1],'FontSize',18,'FontWeight','normal','HorizontalAlignment','Center','Rotation',angle);
        else       
            th = text([loc loc],[pos(3)+nsoffset pos(3)+nsoffset],'n.s.','Color',[0.1 0.1 0.1],'FontSize',16,'FontWeight','normal','HorizontalAlignment','Center','Rotation',angle);
        end
    end
end