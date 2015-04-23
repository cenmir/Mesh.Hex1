classdef Hex1 < matlab.mixin.Copyable
    %Hex1 P1 Hexahedral structured mesh rectangle
    %   o = Hex1(x0,x1,nxe,y0,y1,nye,z0,z1,nze)
    
    properties
        Connectivity
        Points
        Faces
        XC
        YC
        ZC
        edges
        Element
        ElementType
        
        HangNodes
        
        
        nnod
        nele
    end
    
    properties (Hidden)
       
    end
    
    properties (Access = private)
        x0
        x1
        y0
        y1
        z0
        z1
        nxe
        nye
        nze
        nx
        ny
        nz
        
        NeighborsC              % Previously computed 
        
    end
    
    methods
        function o = Hex1(x0,x1,nxe,y0,y1,nye,z0,z1,nze)
            % o = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze)
            o.ElementType = 'Hex1Mesh';
            o.nele = nxe*nye*nze;
            
            if o.nele > 100000
                btn = questdlg(['Warning, ',num2str(o.nele),' elements are about to be created. Do you wish to continue?'], ...
                    '¡Ay! ¡Mucho elementos!!','Si','No, no, nooo', 'No, no, nooo');
                
                switch btn
                    case 'Si'
                        disp('Creating huge amounts of elements...')
                    case 'No, no, nooo'
                        error('No mesh is created!')
                        
                end
            end
            
            o.x0 = x0; o.x1 = x1;
            o.y0 = y0; o.y1 = y1;
            o.z0 = z0; o.z1 = z1;
            
            x = linspace(x0,x1,nxe+1);
            y = linspace(y0,y1,nye+1);
            z = linspace(z0,z1,nze+1);
            [MX, MY, MZ] = meshgrid(x,y,z);
            
            o.Points = [MX(:),MY(:),MZ(:)];
            o.XC = MX(:); o.YC = MY(:); o.ZC = MZ(:);
            
            nx = nxe+1;
            ny = nye+1;
            nz = nze+1;
            
            o.nxe = nxe;
            o.nye = nye;
            o.nze = nze;
            
            o.nele = nxe*nye*nze;
            
            nodes = zeros(o.nele, 8);
            c = 1;
            iel = 1;
            for k = 1:nz
                for j = 1:nx
                    for i = 1:ny
                        
                        if i ~= ny && j ~=nx && k ~= nz
                            nodes(iel,:) = [c,c+1,c+ny+1,c+ny,c+nx*ny,c+nx*ny+ny,c+nx*ny+ny+1,c+nx*ny+1];
                            iel = iel + 1;
                        end
                        c = c+1;
                    end
                end
            end
            
            o.nx = nx;
            o.ny = ny;
            o.nz = nz;
            
            Q = zeros(o.nele*6,4);
            edges1 = zeros(o.nele*12,2);
            
            o.Element(o.nele).faces = [];
            o.Element(o.nele).edges = [];
            
            I = [1     2     4     3     5     8     6     7];
            nodes = nodes(:,I);
            
            lof = 1; loe = 1;
            for iel = 1:o.nele
                upf = lof+5;
                upe = loe+11;
                %Faces
                iface = [nodes(iel,[1,2,4,3]);...
                         nodes(iel,[5,6,8,7]);...
                         nodes(iel,[1,5,6,2]);...
                         nodes(iel,[2,6,8,4]);...
                         nodes(iel,[4,8,7,3]);...
                         nodes(iel,[3,7,5,1])];
                Q(lof:upf,:) = iface;
                o.Element(iel).faces = iface;
                
                %edges
                n = nodes(iel,:);
                
                iedges = [n([1,2]);...
                          n([2,4]);...
                          n([4,3]);...
                          n([3,1]);...
                          n([1,5]);...
                          n([3,7]);...
                          n([4,8]);...
                          n([2,6]);...
                          n([6,5]);...
                          n([5,7]);...
                          n([7,8]);...
                          n([8,6])];
                edges1(loe:upe,:) = iedges;
                o.Element(iel).edges = iedges;
                
                
                %counter
                lof = upf +1;
                loe = upe +1;
                
                
            end
            o.Connectivity = nodes;
            
            o.Faces = Q;
            
            o.nnod = length(o.XC);
            o.edges = edges1;
            
            o.HangNodes = [];
            o.Element(1).HangNodes = [];
        end
        
        
        [neighs,o] = Neighbors(o,varargin);
        
        [HangNodes, HangNodesM, o] = RefineLocal(o, ele);
        
        h = vizMesh(o,varargin);
        
        [fi, fix, fiy, fiz, vol] = baseFun(o,iel,X);
        
        [BoundVerts, BoundEle] = boundaryInds(o, varargin);


    end
    
end

