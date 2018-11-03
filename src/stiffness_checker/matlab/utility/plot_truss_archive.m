
function plot_truss(xn,f,Idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp,rho,rho_min,iplot,filename,save_file_path)

display_factor=0.1;
eps=1e-4;
file_format = 'png';
% 'jpeg'
% 'png'
% 'svg'
% 'epsc' (Encapsulated PostScript (EPS) Level 3 color)
% for more format please check here: https://www.mathworks.com/help/matlab/ref/saveas.html

if nargin < 16
    save_file_path = './ComputationRecords/';
    if ~exist(save_file_path, 'dir')
        save_file_path = './';
    end
end

% characteristic distances
xmax=max(max(xn(1,:)));
xmin=min(min(xn(1,:)));

if (nsd >1)
    ymax=max(max(xn(2,:)));
    ymin=min(min(xn(2,:)));
    
    Lcar=max([xmax-xmin;ymax-ymin]);
else
    Lcar=xmax-xmin;
end

if nsd < 3
    h = figure;
    axis equal;
    title('Undeformed mesh and BCs (nsd < 3)');
    xlabel('x axis');
    ylabel('y axis');
    zlabel('z axis');
    hold on;
    plot_mesh_underformed(nel,ien,xn,nnp,nsd,rho,rho_min,iplot);
    % numbers(nel,ien,xn,nnp,nsd);
    plot_bc_displacements(Lcar,display_factor,nnp,Idb,xn,nsd);
    plot_bc_force(Lcar,f,display_factor,nnp,xn,nsd);
    hold off;
    drawnow;
    saveas(h,strcat(save_file_path,filename(1:end-4)),file_format)
else
    h = figure;
    axis equal;
    title('Undeformed mesh and BCs (3D)');
    xlabel('x axis');
    ylabel('y axis');
    zlabel('z axis');
    hold on;
    plot_mesh_underformed3(nel,ien,xn,nnp,nsd,rho,rho_min,iplot);
    % numbers(nel,ien,xn,nnp,nsd);
    plot_bc_displacements(Lcar,display_factor,nnp,Idb,xn,nsd);
    plot_bc_force(Lcar,f,display_factor,nnp,xn,nsd);
    view(26,10.35);
    %view([1,1,1])
    hold off;
    drawnow;
    saveas(h,strcat(save_file_path,filename(1:end-4)),file_format)
end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% undeformed configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_mesh_underformed(nel,ien,xn,nnp,nsd,rho,rho_min,iplot)

scale = 1;

for e=1:nel
    if rho(e) > rho_min
        node1=ien(1,e);
        node2=ien(2,e);
        x0=[xn(1,node1);xn(1,node2)];
        if (nsd > 1)
            y0=[xn(2,node1);xn(2,node2)];
        else
            y0=[0;0];
        end;
        
        plot(x0,y0,'b-o','LineWidth',rho(e)/scale,'Color',[1-rho(e),0,rho(e)]);
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% undeformed configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_mesh_underformed3(nel,ien,xn,nnp,nsd,rho,rho_min,iplot)

scale = 5;

for e=1:nel
    if rho(e) > rho_min
        node1=ien(1,e);
        node2=ien(2,e);
        x0=[xn(1,node1);xn(1,node2)];
        y0=[xn(2,node1);xn(2,node2)];
        z0=[xn(3,node1);xn(3,node2)];
        plot3(x0,y0,z0,'b','LineWidth',rho(e)/scale,'Color',[1-rho(e),0,rho(e)]);
    end
end
return

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % numbers the nodes and elements - truss and beam %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function numbers(nel,ien,xn,nnp,nsd)
% % number the nodes in the undeformed configuration
% for n=1:nnp
%     if (nsd > 1)
%         text(xn(1,n),xn(2,n), num2str(n),'FontSize',14);
%     else
%         text(xn(1,n),0, num2str(n),'FontSize',14);
%     end;
% end
% % number the elements in the undeformed configuration
% for e=1:nel
%     node1=ien(1,e);
%     node2=ien(2,e);
%     %compute the coordinates of the middle point
%     xg=0.5*(xn(:,node1)+xn(:,node2));
%     s=sprintf('(%d)', e);
%     if (nsd > 1)
%         text(xg(1),xg(2), s,'FontSize',14);
%     else
%         text(xg(1),0, s,'FontSize',14);
%     end;
% end
% end

function plot_bc_displacements(Lcar,display_factor,nnp,Idb,xn,nsd)
    alpha=display_factor*Lcar; % scale factor  for bc symbols
    
    for P=1:nnp
        switch nsd
            case 1
                if (Idb(1,P) ~= 0)
                    bc_symbols([xn(1,P),0],alpha,2);
                end
            case 2
                if ((Idb(1,P) ~= 0) & (Idb(2,P) ~= 0))
                    bc_symbols(xn(:,P),alpha,3);
                end
                if ((Idb(1,P) ~= 0) & (Idb(2,P) == 0))
                    bc_symbols(xn(:,P),alpha,2);
                end
                if ((Idb(1,P) == 0) & (Idb(2,P) ~= 0))
                    bc_symbols(xn(:,P),alpha,1);
                end
            case 3
                % added by Yijiang
                if (Idb(1,P) ~= 0) || (Idb(2,P) ~= 0) || (Idb(3,P) ~= 0)
                    bc_symbols(xn(:,P),alpha,7, Idb(1,P), Idb(2,P), Idb(3,P));
                end
            otherwise
                error('ERROR: dof bigger than three!');
        end
    end
    return
    
    
    function bc_symbols(xp,alpha,symbol, Tx, Ty, Tz)
    
    switch symbol
        case 1
            % v fixed
            x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
            y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3,
                circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(7/8)*alpha],alpha/8);
            end;
            
            
        case 2
            % u fixed
            x=[xp(1);xp(1)-(3/4)*alpha;xp(1)-(3/4)*alpha;xp(1)];
            y=[xp(2);xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha;xp(2)];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3,
                circle([xp(1)-(7/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
            end;
            
        case 3
            % u and v fixed
            x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
            y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3,
                line([xp(1)-(alpha/4)+i*alpha/4;xp(1)-(alpha/2)+i*(alpha/4)], ...
                    [xp(2)-(3/4)*alpha;xp(2)-alpha],'Color','k','LineWidth',1.2);
            end;
            
        case 4
            % v and theta fixed
            x=[xp(1)-alpha/2;xp(1)+alpha/2];
            y=[xp(2);xp(2)];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3,
                circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(1/8)*alpha],alpha/8);
            end;
            
        case 5
            % u and theta fixed
            x=[xp(1);xp(1)];
            y=[xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3,
                circle([xp(1)-(1/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
            end;
            
        case 6
            % u, v and theta fixed
            line([xp(1)-alpha/2;xp(1)+alpha/2],[xp(2),xp(2)],'Color','k','LineWidth',1.2);
            for i=0:3,
                line([xp(1)-alpha/2+(i+1)*alpha/4, xp(1)-alpha/2+i*alpha/4],[xp(2),xp(2)-alpha/4]...
                    ,'Color','k','LineWidth',1.2);
            end;
        case 7
            % 3d case, added by Yijiang
            % u, v, and w fixed
            x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
            y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
%             z=[xp(3);xp(3)-(3/4)*alpha;xp(3)-(3/4)*alpha;xp(3)];
            
%             line(x,y,'Color','k','LineWidth',1.2);
            scatter3(xp(1),xp(2),xp(3),'filled', 'k');

            if Tx
                plot3([xp(1)-(alpha); xp(1)],...
                    [xp(2); xp(2)],...
                    [xp(3); xp(3)],'Color','k','LineWidth',1.2);
            end
            
            if Ty
                plot3([xp(1); xp(1)],...
                    [xp(2)-(alpha); xp(2)],...
                    [xp(3); xp(3)],'Color','k','LineWidth',1.2);
            end
            
            if Tz
                plot3([xp(1); xp(1)],...
                    [xp(2); xp(2)],...
                    [xp(3)-(alpha); xp(3)],'Color','k','LineWidth',1.2);
            end

    end
    return
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % boundary conditions on force %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_bc_force(Lcar,f,display_factor,nnp,xn,nsd)
            delta=display_factor*Lcar/max(max(abs(f))); % scale factor for force b.c.
            alpha=2*display_factor*Lcar; % scale factor for moment
            
            for N=1:nnp
                
                if (nsd > 2)
                    if ( (f(1,N) ~=0 ) | (f(2,N) ~= 0)  | (f(3,N) ~= 0))
                        quiver3(xn(1,N),xn(2,N),xn(3,N),f(1,N),f(2,N),f(3,N),delta,'k','LineWidth',1.5);
                    end;
                elseif ( nsd >1)
                    if ( (f(1,N) ~=0 ) | (f(2,N) ~= 0))
                        quiver(xn(1,N),xn(2,N),f(1,N),f(2,N),delta,'k','LineWidth',1.5);
                    end;
                else
                    if (f(1,N) ~=0) quiver(xn(1,N),0,f(1,N),0,delta,'k');
                    end;
                end;
            end;
            return
            
            function circle(x0,r);
                theta=0:0.1:2*pi;
                x=r*cos(theta)+x0(1);
                y=r*sin(theta)+x0(2);
                
                plot(x,y,'k','LineWidth',1.2);
                return