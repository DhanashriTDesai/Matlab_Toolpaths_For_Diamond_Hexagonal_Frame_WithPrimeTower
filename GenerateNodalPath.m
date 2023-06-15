%% GenerateNodalPath function
%%% This function generates a tool-path, according to no. of cells and nodal data,
% by defining a nodal-path such that entire 2D-frame structure is traced with minimal no. of brakes %%%

function [NodalPath] = GenerateNodalPath(CellType,nnodes,nCellH,nCellV)

NodalPath=[]; k=1;
switch CellType
    case 'StretchDominatedDiamond'                 
                    RectangularDomainNodalPath=[nnodes+1 nnodes+2 nnodes+3 nnodes+4 nnodes+1];
                    for i=1:size(RectangularDomainNodalPath,2)
                        NodalPath(k)=RectangularDomainNodalPath(i); k=k+1;
                    end
                    
                    HorizontalNodalPath=zeros(nCellH-1,nCellV); %NumOfHorizontalNodalPaths=nCellH-1 and LengthOfHorizontalPath=nCellV
                    for i=1:(nCellH-1)
                        for j=1:nCellV
                            HorizontalNodalPath(i,j)=2*i*nCellV+(j+i);
                            NodalPath(k)=HorizontalNodalPath(i,j); k=k+1;
                        end
                    end
                    
                    %  generalized zigzag path 
                    for itr=0:(nCellH-1)
                        if(rem(itr,2)==0)
                            i = itr*(2*nCellV+1)+1;
                            j = i+nCellV;
                            b = j+nCellV;
                            a = j-1;
                            while(a>=i || b>=j)
                                if(b>=j)
                                    NodalPath(k)=b; k=k+1;
                                end
                                if(a>=i)
                                    NodalPath(k)=a; k=k+1;
                                end
                                a=a-1; 
		                        b=b-1;
                            end
                            b = b+2;
                            while(b<=j+nCellV)
                                NodalPath(k)=b; k=k+1;
                                b=b+1;
                            end
                            i = j;
                            j = i+nCellV+1;
                            a = j-2;
                            b = j+nCellV-1;
                            while(a>=i || b>=j)
                                if(b>=j)
                                    NodalPath(k)=b; k=k+1;
                                end
                                if(a>=i)
                                    NodalPath(k)=a; k=k+1;
                                end
                                a=a-1;
                                b=b-1;
                            end
                        else
                            i = itr*(2*nCellV+1)+1;
                            j = i+nCellV;
                            a = i;
                            b = j;
                            while(a<i+nCellV || b<j+nCellV)
                                if(b<j+nCellV)
                                    NodalPath(k)=b; k=k+1;
                                end
                                if(a<i+nCellV)
                                    NodalPath(k)=a; k=k+1;
                                end
                                a=a+1;
                                b=b+1;
                            end
                            b = b;
                            while(j<b)
                                NodalPath(k)=b; k=k+1;
                                b=b-1;
                            end
                            i = j;
                            j = i+nCellV+1;
                            a = i;
                            b = j;
                            while(a<i+nCellV+1 || b<j+nCellV)
                                if(a<i+nCellV+1)
                                    NodalPath(k)=a; k=k+1;
                                end
                                if(b<j+nCellV)
                                    NodalPath(k)=b; k=k+1;
                                end
                                a=a+1;
                                b=b+1;
                            end
                        end
                    end
                    
    case 'BendingDominatedHexagon'
                    if (mod(nCellH,2)==0)
                        RectangularDomainNodalPath=[1 (nnodes+2) (nnodes+1) (nnodes-nCellH) 1];
                    else
                        RectangularDomainNodalPath=[1 (nCellH+1) nnodes (nnodes-nCellH) 1];
                    end
                    
                    for i=1:size(RectangularDomainNodalPath,2)
                        NodalPath(k)=RectangularDomainNodalPath(i); k=k+1;
                    end
                    
                    % generalized zigzag path 
                    if (mod(nCellH,2)==0)
                            i = 0; j = 0; jump = fix(nCellH/2)*2+1;
                            for itr=1:nCellV
                                i = i+jump;
                                a = i+1;
                                b = i-jump+2;
                                alim = i+jump-1;
                                blim = i;
                        
                                while(a<=alim || b<=blim)
                                    if(a<=alim)
                                        NodalPath(k) = a; k=k+1;
                                    end
                                    a=a+1;
                                    if(a<=alim)
                                        NodalPath(k) = a; k=k+1;
                                    end
                                   
                                    if(b<=blim)
                                        NodalPath(k) = b; k=k+1;
                                    end
                                    b=b+1;
                                    if(b<=blim)
                                        NodalPath(k) = b; k=k+1;
                                    end
                                    a=a+1;
                                    b=b+1;
                                end
                        
                                i = i+jump;
                                a = i+jump;
                                b = i;
                                alim = i+1;
                                blim = i-jump+1;
                                
                                while(a>=alim || b>=blim)
                                    if(b>=blim)
                                        NodalPath(k) = b; k=k+1;
                                    end
                                    b=b-1;
                                    if(a>=alim)
                                        NodalPath(k) = a; k=k+1;
                                    end
                                    a=a-1;
                                    if(a>=alim)
                                        NodalPath(k) = a; k=k+1;
                                    end    
                                    if(b>=blim)
                                        NodalPath(k) = b; k=k+1;
                                    end
                                    a=a-1;
                                    b=b-1;
                                end
                            end

                    else
                                i = 0; j = 0; jump = fix(nCellH/2)*2+2;
                                for itr=1:nCellV
                                    i = i+jump;
                                    a = i+1;
                                    b = i-jump+2;
                            
                                    while(a<=i+jump && b<=i)
                                        NodalPath(k) = a; k=k+1;
                                        a=a+1;
                                        if(a<=i+jump)
                                            NodalPath(k) = a; k=k+1;
                                        end
                                        NodalPath(k) = b; k=k+1;
                                        b=b+1;
                                        if(b<=i)
                                            NodalPath(k) = b; k=k+1;
                                        end
                                        a=a+1;
                                        b=b+1;
                                    end

                                    i = i+jump;
                                    a = i+jump;
                                    b = i;
                                    NodalPath(k) = a; k=k+1;
                                    a=a-1;
                                    
                                    while(a>=i+1 && b>=i-jump+1)
                                        NodalPath(k) = b; k=k+1;
                                        b=b-1;
                                        if(b>=i-jump+1)
                                            NodalPath(k) = b; k=k+1;
                                        end
                                        NodalPath(k) = a; k=k+1;
                                        a=a-1;
                                        if(a>=i+1)
                                            NodalPath(k) = a; k=k+1;
                                        end
                                        a=a-1;
                                        b=b-1;
                                    end

                                end
                    end

end
NodalPath=NodalPath';
end