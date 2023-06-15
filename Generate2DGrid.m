%% Generate2DGrid function
%%% This function generates 2D grid containing nodes and elements according to cell parameters and no. of cells %%%

function [LCell,Lh,Lv,nelem,nnodes,nodeID,nx,ny,Nrows,Ncolumns,LhDomain,LvDomain] = Generate2DGrid(CellType,nCellH,nCellV,tCell,rho_rel,theta)

switch CellType
    case 'StretchDominatedDiamond'
                    Nrows = 2*nCellH+1; Ncolumns = 2*nCellV+1; % no of rows and columns in the grid
                    LCell=round(2*sqrt(3)*tCell/rho_rel); % LCell=2*sqrt(3)*tCell/rho_rel theoretical relative density of diamond cell
                    Lh=LCell*cosd(theta); Lv=LCell*sind(theta); % horizontal and vertical components of LCell
                    x0=40+2*nCellV*Lv; y0=60; % start point for the printing 
                    nelem = 0; nnodes = 0;
                    
                    % Co-ordinates of nodes in odd rows
                    for i = 1:2:Nrows 
                        for j = 2:2:Ncolumns
                            nnodes = nnodes + 1;
                            nodeID(i,j) = ((i-1)*Ncolumns + j)/2;
                            nx(nodeID(i,j)) = x0 - Lv * (i-1);
                            ny(nodeID(i,j)) = y0 + Lh * (j-1);
                        end
                    end
                    
                    % Co-ordinates of nodes in even rows
                    for i = 2:2:Nrows 
                        for j = 1:2:Ncolumns
                            nnodes = nnodes + 1;
                            nodeID(i,j) = ((i-1)*Ncolumns + j)/2;
                            nx(nodeID(i,j)) = x0 - Lv * (i-1);
                            ny(nodeID(i,j)) = y0 + Lh * (j-1);        
                        end
                    end
                    
                    % Connect horizontal elements in odd rows
                    for i = 1:2:Nrows
                        for j = 2:2:(Ncolumns-2)
                            nelem = nelem + 1;
                            elemcon(nelem,1) = nodeID(i,j);
                            elemcon(nelem,2) = nodeID(i,j+2);
                        end
                    end
                    
                    % Connect horizontal elements in even rows
                    for i = 2:2:Nrows
                        for j = 1:2:(Ncolumns-2)
                            nelem = nelem + 1;
                            elemcon(nelem,1) = nodeID(i,j);
                            elemcon(nelem,2) = nodeID(i,j+2);
                        end
                    end
                            
                    % Connect inclined elements
                    for i = 1:2:(Nrows-1)
                        for j = 1:2:(Ncolumns-2)
                            nelem = nelem + 1;
                            elemcon(nelem,1) = nodeID(i,j+1);
                            elemcon(nelem,2) = nodeID(i+1,j);
                    
                            nelem = nelem + 1;
                            elemcon(nelem,1) = nodeID(i,j+1);
                            elemcon(nelem,2) = nodeID(i+1,j+2);
                        end
                    end       
                    
                    % Connect inclined elements
                    for i = 3:2:Nrows
                        for j = 1:2:(Ncolumns-2)
                            nelem = nelem + 1;
                            elemcon(nelem,1) = nodeID(i,j+1);
                            elemcon(nelem,2) = nodeID(i-1,j);
                    
                            nelem = nelem + 1;
                            elemcon(nelem,1) = nodeID(i,j+1);
                            elemcon(nelem,2) = nodeID(i-1,j+2);  
                        end
                    end
                           
                    % Connect remaining inclined elements
                    if (mod(Nrows,2)==0)
                        for i = 1:2:(Nrows-1)
                            nelem = nelem + 1;
                            elemcon(nelem,1) = nodeID(i,Ncolumns);
                            elemcon(nelem,2) = nodeID(i+1,Ncolumns-1);
                        end
                    
                        for i = 2:2:(Nrows-1)
                            nelem = nelem + 1;  
                            elemcon(nelem,1) = nodeID(i,Ncolumns-1);
                            elemcon(nelem,2) = nodeID(i+1,Ncolumns);  
                        end        
                    end
                    
                    figure(1)
                    for i = 1:nelem
                        node1 = elemcon(i,1); node2 = elemcon(i,2);
                        L(i) = sqrt((nx(node1)-nx(node2))^2 + (ny(node1)-ny(node2))^2);
                        plot([nx(node1) nx(node2)], [ny(node1) ny(node2)],'-r');
                        hold on
                    end
                    
                    % vertices of the rectangular domain
                    nx(nnodes+1)=nx(nnodes); ny(nnodes+1)=ny(1)-Lh; 
                    nx(nnodes+2)=nx(1); ny(nnodes+2)=ny(1)-Lh; 
                    nx(nnodes+3)=nx(1); ny(nnodes+3)=ny(nnodes)+Lh;
                    nx(nnodes+4)=nx(nnodes); ny(nnodes+4)=ny(nnodes)+Lh;
                    LhDomain=nx(nnodes+2)-nx(nnodes+1); LvDomain=ny(nnodes+4)-ny(nnodes+1);
                    
                    for i = 1:length(nx)
                        text(nx(i),ny(i),num2str(i),'Fontsize',14,'Color','k');
%                         text(nx(i),(ny(i)-5),num2str([nx(i) ny(i)]),'Fontsize',10,'Color','b');
                        plot(nx(i),ny(i),'b*',MarkerSize=4);
                        plot(nx(i),ny(i),'bo',MarkerSize=4);
                    end
                    axis equal; xlim([0 230]); ylim([0 230]); 

    case 'BendingDominatedHexagon'
            Nrows = 2*nCellV+1; Ncolumns = 2*nCellH+2; % no of rows and columns in the grid
            LCell=round(2/sqrt(3)*tCell/rho_rel);
            Lh = LCell*cosd(theta); Lv = LCell*sind(theta); % horizomtal and vertical lengths
            x0=40; y0=60; % start point for the printing 
            nelem = 0; nnodes = 0;
            
            % Co-ordinates of nodes in odd rows
            for i = 1:2:Nrows 
                for j = 4:4:Ncolumns
                    nnodes = nnodes + 1;
                    nodeID(i,j) = ((i-1)*Ncolumns+j)/2;
                    nx(nodeID(i,j)) = x0+(j/2)*Lh + (j/2-1)*LCell;
                    ny(nodeID(i,j)) = y0+Lv * (i-1);
                end
            
                for j = 5:4:Ncolumns
                    nnodes = nnodes + 1;
                    nodeID(i,j) = ((i-1)*Ncolumns+j+1)/2;
                    nx(nodeID(i,j)) = x0+(j-1)/2*(Lh+LCell);
                    ny(nodeID(i,j)) = y0+Lv * (i-1);
                end
            
                nnodes = nnodes + 1;
                nodeID(i,1) = (i-1)/2*Ncolumns + 1;
                nx(nodeID(i,1)) = x0;
                ny(nodeID(i,1)) = y0+Lv * (i-1);
            end
            
            % Co-ordinates of nodes in even rows
            for i = 2:2:Nrows 
                for j = 2:4:Ncolumns
                    nnodes = nnodes + 1;
                    nodeID(i,j) = ((i-1)*Ncolumns+j)/2;
                    nx(nodeID(i,j)) = x0+(j/2)*Lh + (j-2)/2*LCell;
                    ny(nodeID(i,j)) = y0+Lv * (i-1);        
                end
            
                for j = 3:4:Ncolumns
                    nnodes = nnodes + 1;
                    nodeID(i,j) = ((i-1)*Ncolumns+j+1)/2;
                    nx(nodeID(i,j)) = x0+(j-1)/2*(Lh+LCell);
                    ny(nodeID(i,j)) = y0+Lv * (i-1);        
                end
            end
            
            % additional nodes for tracing the domain
            if (mod(nCellH,2)==0)
                nx(nnodes+1)=nx(nnodes)+Lh; ny(nnodes+1)=ny(nnodes);
                nx(nnodes+2)=nx(nCellH+1)+Lh; ny(nnodes+2)=ny(1); 
                LhDomain=nx(nnodes+2)-nx(1); LvDomain=ny(nnodes+1)-ny(1); k=nnodes+2;
            else
                LhDomain=nx(nCellH+1)-nx(1); LvDomain=ny(nnodes)-ny(1); k=nnodes;
            end
            
            % Connect horizontal elements in odd rows
            for i = 1:2:Nrows
                for j = 4:4:(Ncolumns-1)
                    nelem = nelem + 1;
                    elemcon(nelem,1) = nodeID(i,j);
                    elemcon(nelem,2) = nodeID(i,j+1);
                end
            end
            
            % Connect horizontal elements in even rows
            for i = 2:2:Nrows
                for j = 2:4:(Ncolumns-1)
                    nelem = nelem + 1;
                    elemcon(nelem,1) = nodeID(i,j);
                    elemcon(nelem,2) = nodeID(i,j+1);
                end
            end
                    
            % Connect inclined elements in odd rows
            for i = 1:2:(Nrows-1)
                for j = 1:4:(Ncolumns-1)
                    nelem = nelem + 1;
                    elemcon(nelem,1) = nodeID(i,j);
                    elemcon(nelem,2) = nodeID(i+1,j+1);
                end
            
                for j = 4:4:Ncolumns
                    nelem = nelem + 1;
                    elemcon(nelem,1) = nodeID(i,j);
                    elemcon(nelem,2) = nodeID(i+1,j-1);
                end
            end  
            
            % Connect inclined elements in even rows
            for i = 2:2:Nrows
                for j = 2:4:Ncolumns
                    nelem = nelem + 1;
                    elemcon(nelem,1) = nodeID(i,j);
                    elemcon(nelem,2) = nodeID(i+1,j-1);
                end
            
                for j = 3:4:(Ncolumns-1)
                    nelem = nelem + 1;
                    elemcon(nelem,1) = nodeID(i,j);
                    elemcon(nelem,2) = nodeID(i+1,j+1);  
                end
            end
            
            % Plot the ground structure - elements + nodes
            figure(1)
            for i = 1:nelem
                node1 = elemcon(i,1); node2 = elemcon(i,2);
                L(i) = sqrt((nx(node1)-nx(node2))^2 + (ny(node1)-ny(node2))^2);
                plot([nx(node1) nx(node2)], [ny(node1) ny(node2)],'-r');
                hold on
            end
            
            for i = 1:k
                text(nx(i),ny(i),num2str(i),'Fontsize',14,'Color','k');
%                 text(nx(i)-15,ny(i)-5,num2str([nx(i) ny(i)]),'Fontsize',10,'Color','b');
                plot(nx(i),ny(i),'b*',MarkerSize=4);
                plot(nx(i),ny(i),'bo',MarkerSize=4);
            end
            axis equal; xlim([0 270]); ylim([0 270]);
end
end