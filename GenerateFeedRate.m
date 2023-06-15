%% GenerateFeedRate function
%%% This function generates the feed rate data according to cell parameters, no. of cells and no. of layers
% by incrementing the extrusion co-rdinate (E) through constant value %%%

function [E,nExtrusion] = GenerateFeedRate(CellType,E0,nelem,dLdE,LCell,Lv,LhDomain,LvDomain,nCellH,nCellV,Nrows,nLayer)

switch CellType
    case 'StretchDominatedDiamond'
                % Generate feed rate data for print path
                nExtrusion=nelem-2*(nCellV-1)+4; % no of extrusion moves
                
                for i=1:nLayer
                    for j=1:nExtrusion
                        if i==1 && j==1
                            E((i-1)*nExtrusion+j)=E0+LhDomain/dLdE;
                        elseif j==1 || j==3
                            E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+LhDomain/dLdE;
                        elseif j==2 || j==4
                            E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+LvDomain/dLdE;
                        else 
                            E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+LCell/dLdE;
                        end
                    end
                end

    case 'BendingDominatedHexagon'
                % Generate feed rate data for print path
                nHorizontalDoubleExtrusion=0;
                if mod(nCellH,2)==0
                    for i=2:(Nrows-1)
                        nHorizontalDoubleExtrusion=nHorizontalDoubleExtrusion+fix(nCellH/2);
                    end
                else
                    for i=2:(Nrows-1)
                        if mod(i,2)==0
                            nHorizontalDoubleExtrusion=nHorizontalDoubleExtrusion + (fix(nCellH/2)+1);
                        else
                            nHorizontalDoubleExtrusion=nHorizontalDoubleExtrusion + fix(nCellH/2);
                        end
                    end
                end
                
                if(mod(nCellH,2)==1)
                    nExtrusion=nelem+4+nHorizontalDoubleExtrusion+nCellV; % no of extrusion moves
                    jump(1)=5+(2*nCellH+1);
                    for i=2:nCellV
                        jump(i) = jump(i-1) + 2*(2*nCellH+1) + 1;
                    end
                    for i=1:nLayer
                        for j=1:nExtrusion
                            if i==1 && j==1
                                E((i-1)*nExtrusion+j)=E0+LhDomain/dLdE;
                            elseif j==1 || j==3
                                E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+LhDomain/dLdE;
                            elseif j==2 || j==4
                                E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+LvDomain/dLdE;
                            elseif ismember(j,jump)
                                E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+ (2*Lv/dLdE);
                            else 
                                E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+LCell/dLdE;
                            end
                        end
                    end

                else
                    nExtrusion=nelem+4+nHorizontalDoubleExtrusion; % no of extrusion moves
                    for i=1:nLayer
                        for j=1:nExtrusion
                            if i==1 && j==1
                                E((i-1)*nExtrusion+j)=E0+LhDomain/dLdE;
                            elseif j==1 || j==3
                                E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+LhDomain/dLdE;
                            elseif j==2 || j==4
                                E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+LvDomain/dLdE;
                            else 
                                E((i-1)*nExtrusion+j)=E((i-1)*nExtrusion+j-1)+LCell/dLdE;
                            end
                        end
                    end
                end

end
E=E';
end