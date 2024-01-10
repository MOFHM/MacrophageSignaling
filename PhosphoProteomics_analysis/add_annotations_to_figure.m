function[] = add_annotations_to_figure(fig,x_data,y_data,annotation_vector,nr_annotations)

%Parameters:
dim_frame = 20000;
point_size = round(0.0095*dim_frame);
dim_text_x = round(0.015*dim_frame);
dim_text_y = round(0.045*dim_frame);

% point_size = round(0.0055*dim_frame);
% dim_text_x = round(0.010*dim_frame);
% dim_text_y = round(0.025*dim_frame);

ax_lim_x = fig.CurrentAxes.XLim;
ax_lim_y = fig.CurrentAxes.YLim;
theta = 0 : 0.1 : 2*pi;
radius = [0.2,0.25,0.17,0.15,0.13,0.1,0.07];

annotation_score = sqrt(2.*(y_data).^2+x_data'.^2);
annotation_score =  annotation_score.*sign(x_data)';
[a1,b1] = sort(annotation_score,'descend');
a11 = isfinite(a1);
ind1 = find(a11,1,'first');
if(length(b1)<nr_annotations/2)
    nr_annotations = round(length(b1)/3);
end
frame = ones(dim_frame,dim_frame);
[normx, normy] = coord2norm(gca, [x_data],[y_data]);
normx = round(normx.*dim_frame);
ic = find(isinf(normx)==0);
normx = normx(ic);
normy = normy(ic);
normy = dim_frame - round(normy.*dim_frame);
ic = find(isinf(normy)==0);
normx = normx(ic);
normy = normy(ic);
ic = find(isnan(normx));
normx(ic) = [];
normy(ic) = [];
ic = find(isnan(normy));
normx(ic) = [];
normy(ic) = [];
for rr = 1:length(normx)
    c1 = normx(rr)>2*point_size && normy(rr)>2*point_size;
    c2 = normx(rr)<dim_frame-2*point_size && normy(rr)<dim_frame-2*point_size;
    if(c1 && c2)
        frame(normy(rr)-point_size:normy(rr)+point_size,normx(rr)-point_size:normx(rr)+point_size) = 0;
    else if(c1==1 && c2==0)
            frame(normy(rr)-point_size:normy(rr),normx(rr)-point_size:normx(rr)) = 0;
        else if(c1==0 && c2==1)
                frame(normy(rr):normy(rr)+point_size,normx(rr):normx(rr)+point_size) = 0;
            end
        end
    end
end

for an = 0:nr_annotations-1
    x_b = x_data(b1(ind1+an))-0.005;
    x_e = x_data(b1(ind1+an))-radius(1) * sin(theta(round(length(theta)/4)));
    y_b = y_data(b1(ind1+an))+0.005;
    y_e = y_data(b1(ind1+an))+radius(1) * cos(theta(round(length(theta)/4)));
    if( (y_e>0.98*max(ax_lim_y) || y_e<0.98*min(ax_lim_y)) || (x_e>0.98*max(ax_lim_x) || x_e<0.98*min(ax_lim_x)) )
        x_b = 0.5;
        x_e = 0.5;
        y_b = 1;
        y_e = 1;
    end
    [normx, normy] = coord2norm(gca, [x_e,x_b],[y_e,y_b]);
    normx = round(normx(1)*dim_frame);
    normy = dim_frame-round(normy(1)*dim_frame);
    overlap_value =[];
    for ixj = 1:length(radius)
        nr_it = 1;
        while(nr_it<length(theta)-1)
            % while (length(find(frame(normy-dim_text_x:normy+dim_text_x,normx-dim_text_y:normx+dim_text_y)))<0.85*4*(dim_text_y*dim_text_x))
            %if (nr_it==length(theta)) break; end
            nr_it = nr_it+1;
            x_b = x_data(b1(ind1+an))-0.005;
            x_e = x_data(b1(ind1+an))-radius(ixj) * sin(theta(nr_it));
            y_b = y_data(b1(ind1+an))+0.005;
            y_e = y_data(b1(ind1+an))+radius(ixj) * cos(theta(nr_it));
            if( (y_e>0.98*max(ax_lim_y) || y_e<0.98*min(ax_lim_y)) || (x_e>0.98*max(ax_lim_x) || x_e<0.98*min(ax_lim_x)) ) continue; end
            [normx, normy] = coord2norm(gca, [x_e,x_b],[y_e,y_b]);
            normx = round(normx(1)*dim_frame);
            normy = dim_frame-round(normy(1)*dim_frame);
            if(normy>dim_text_x && normx>dim_text_y && normy+dim_text_x < size(frame,2) && normx+dim_text_y < size(frame,1))
                overlap_value = [overlap_value;[length(find(frame(normy-dim_text_x:normy+dim_text_x,normx-dim_text_y:normx+dim_text_y))),ixj,nr_it]];
            end
        end
    end
    
    ic1 = find(overlap_value(:,1)==max(overlap_value(:,1)));
    ic1 = ic1(1);
    x_b = x_data(b1(ind1+an))-0.005;
    x_e = x_data(b1(ind1+an))-radius(overlap_value(ic1,2)) * sin(theta(overlap_value(ic1,3)));
    y_b = y_data(b1(ind1+an))+0.005;
    y_e = y_data(b1(ind1+an))+radius(overlap_value(ic1,2)) * cos(theta(overlap_value(ic1,3)));
    [normx, normy] = coord2norm(gca, [x_e,x_b],[y_e,y_b]);
    normx = round(normx(1)*dim_frame);
    normy = dim_frame-round(normy(1)*dim_frame);
    
    if( (y_e>0.98*max(ax_lim_y) || y_e<0.98*min(ax_lim_y)) || (x_e>0.98*max(ax_lim_x) || x_e<0.98*min(ax_lim_x)) ) continue; end
    if (length(find(frame(normy-dim_text_x:normy+dim_text_x,normx-dim_text_y:normx+dim_text_y))) < 0.85*4*(dim_text_y*dim_text_x))
        continue
    end
    
    Annotate(gca,'textarrow',[x_e,x_b],[y_e,y_b],'String',annotation_vector(b1(ind1+an)),'FontSize',10,'HeadLength',9,'TextBackgroundColor','none','TextEdgeColor','none');
    [normx, normy] = coord2norm(gca, [x_e,x_b],[y_e,y_b]);
    normx = round(normx(1)*dim_frame);
    normy = dim_frame-round(normy(1)*dim_frame);
    if(normx>2*dim_text_y && normy>2*dim_text_x)
        frame(normy-dim_text_x:normy+dim_text_x,normx-dim_text_y:normx+dim_text_y) = 0;
    else
        frame(normy:normy+dim_text_x,normx:normx+dim_text_y) = 0;
    end
end

ind1 = find(a11,1,'last');
for an = 0:nr_annotations-1
    x_b = x_data(b1(ind1-an))-0.005;
    x_e = x_data(b1(ind1-an))-radius(1) * sin(theta(round(length(theta)/4)));
    y_b = y_data(b1(ind1-an))+0.005;
    y_e = y_data(b1(ind1-an))+radius(1) * cos(theta(round(length(theta)/4)));
    if( (y_e>0.98*max(ax_lim_y) || y_e<0.98*min(ax_lim_y)) || (x_e>0.98*max(ax_lim_x) || x_e<0.98*min(ax_lim_x)) )
        x_b = 0.5;
        x_e = 0.5;
        y_b = 1;
        y_e = 1;
    end
    [normx, normy] = coord2norm(gca, [x_e,x_b],[y_e,y_b]);
    normx = round(normx(1)*dim_frame);
    normy = dim_frame-round(normy(1)*dim_frame);
    overlap_value =[];
    for ixj = 1:length(radius)
        nr_it = 1;
        while(nr_it<length(theta)-1)
            %while (length(find(frame(normy-dim_text_x:normy+dim_text_x,normx-dim_text_y:normx+dim_text_y)))<0.85*4*(dim_text_y*dim_text_x))
            nr_it = nr_it+1;
            %if (nr_it==length(theta)) break; end
            x_b = x_data(b1(ind1-an))-0.005;
            x_e = x_data(b1(ind1-an))-radius(ixj) * sin(theta(nr_it));
            y_b = y_data(b1(ind1-an))+0.005;
            y_e = y_data(b1(ind1-an))+radius(ixj) * cos(theta(nr_it));
            if( (y_e>0.98*max(ax_lim_y) || y_e<0.98*min(ax_lim_y)) || (x_e>0.98*max(ax_lim_x) || x_e<0.98*min(ax_lim_x)) ) continue; end
            [normx, normy] = coord2norm(gca, [x_e,x_b],[y_e,y_b]);
            normx = round(normx(1)*dim_frame);
            normy = dim_frame-round(normy(1)*dim_frame);
            if(normy>dim_text_x && normx>dim_text_y && normy+dim_text_x < size(frame,2) && normx+dim_text_y < size(frame,1))
                overlap_value = [overlap_value;[length(find(frame(normy-dim_text_x:normy+dim_text_x,normx-dim_text_y:normx+dim_text_y))),ixj,nr_it]];
            end
        end
    end
    
    ic1 = find(overlap_value(:,1)==max(overlap_value(:,1)));
    ic1 = ic1(1);
    x_b = x_data(b1(ind1-an))-0.005;
    x_e = x_data(b1(ind1-an))-radius(overlap_value(ic1,2)) * sin(theta(overlap_value(ic1,3)));
    y_b = y_data(b1(ind1-an))+0.005;
    y_e = y_data(b1(ind1-an))+radius(overlap_value(ic1,2)) * cos(theta(overlap_value(ic1,3)));
    [normx, normy] = coord2norm(gca, [x_e,x_b],[y_e,y_b]);
    normx = round(normx(1)*dim_frame);
    normy = dim_frame-round(normy(1)*dim_frame);
    
    if( (y_e>0.98*max(ax_lim_y) || y_e<0.98*min(ax_lim_y)) || (x_e>0.98*max(ax_lim_x) || x_e<0.98*min(ax_lim_x)) ) continue; end
    if (length(find(frame(normy-dim_text_x:normy+dim_text_x,normx-dim_text_y:normx+dim_text_y)))) <0.85*4*(dim_text_y*dim_text_x)
        continue
    end
    Annotate(gca,'textarrow',[x_e,x_b],[y_e,y_b],'String',annotation_vector(b1(ind1-an)),'FontSize',10,'HeadLength',9,'TextBackgroundColor','none','TextEdgeColor','none');
    [normx, normy] = coord2norm(gca, [x_e,x_b],[y_e,y_b]);
    normx = round(normx(1)*dim_frame);
    normy = dim_frame-round(normy(1)*dim_frame);
    if(normx>2*dim_text_y && normy>2*dim_text_x)
        frame(normy-dim_text_x:normy+dim_text_x,normx-dim_text_y:normx+dim_text_y) = 0;
    else
        frame(normy:normy+dim_text_x,normx:normx+dim_text_y) = 0;
    end
end

end