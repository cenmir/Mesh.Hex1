function [HangNodes, HangNodesM, NewElements] = RefineLocal(o, ele)
% RefineLocal Refines local elements by splitting into 8 new elements
% [HangNodes, HangNodesM] = RefineLocal(o, ele)
% [HangNodes, HangNodesM] = o.RefineLocal(ele)
% ele is a list of elements
% list of indices | 'all'
%
% 




if strcmpi(ele,'all')
    ele = 1:size(o.Connectivity,1);
end

%% Refine Hex Mesh Locally
% nodes = o.Connectivity;
% X = o.X;
xnod = o.XC;
ynod = o.YC;
znod = o.ZC;
% Number of nodes, keeps track of the latest node number
nmax = length(xnod);

nodes = o.Connectivity;
X = o.Points;

oldNodes = nodes;

% Loop over all elements that are about to be refined
% For every element we will get 8 subelements
ic = 1;
for iel = ele
    iv = o.Connectivity(iel,:); % local node numbers
    xc = xnod(iv); % local coordinates
    yc = ynod(iv);
    zc = znod(iv);
    
    xm = mean(xc);ym = mean(yc);zm = mean(zc); %mid points of element
    
    edges = o.Element(iel).edges; % Edge list of current elemeno.
    faces = o.Element(iel).faces; % Face list of current elemeno.
    %% Create New Coords
    % We start by creating a list of new points. The list includes points on
    % the mid points of edges, faces and elemeno. 19 points total
    
    % Edge indices, the order is chosen to create the new points in the same
    % order as the mother nodes are numbered (along positive y, pos x and
    % then pos z. See http://www.mirzacenanovic.com/wp-content/uploads/2015/01/2015-01-09-14_50_26-Figure-1_-Hex1Mesh.png
    eind = [1,4,2,3,5,8,6,7,9,10,12,11];
    faceind = [1,3,6,4,5,2];
    %The order in which the subelement points are created is the same order
    %as for the mother element
    XNeind = [1,2,4,5,6,8,12,14,15,16,18,19]; %edge inds
    XNMind = 10; %mid
    XNfind = [3,7,9,11,13,17]; %face
    
    % Notice that we're creating 19 points regardless if these points exists
    % in neighboring elements. We'll deal with that later.
    XN = zeros(19,3);   % 19 New points
    XN(XNeind,:) = (X(edges(eind,1),:)+X(edges(eind,2),:))/2;
    XN(XNMind,:) = [xm,ym,zm];
    XN(XNfind,:)=(X(faces(faceind,1),:)+X(faces(faceind,2),:)+X(faces(faceind,3),:)+X(faces(faceind,4),:))/4;
    
    %Local node numbering same as for mother elements
    locnodes = [1 2 4 5 10 11 13 14;...
                2 3 5 6 11 12 14 15;...        
                4 5 7 8 13 14 16 17;...
                5 6 8 9 14 15 17 18;...
                10 11 13 14 19 20 22 23;...
                11 12 14 15 20 21 23 24;...
                13 14 16 17 22 23 25 26;...
                14 15 17 18 23 24 26 27];
    
    %% Viz
%     xfigure(43435)
%     fele = [6*iel-5;6*iel-4;6*iel-3;6*iel-2;6*iel-1;6*iel-0;];
%     h.patch = patch(o.XC(o.Faces(fele(:),:)'),o.YC(o.Faces(fele(:),:)'),o.ZC(o.Faces(fele(:),:)'),'w','FaceColor','none');
%     xlabel('X'); ylabel('Y'); zlabel('Z')
%     axis equal tight
%     for i = iv
%        text(xnod(i),ynod(i),znod(i),num2str(i),'BackgroundColor','w') 
%     end
%     
%     for i = 1:length(XN)
%         text(XN(i,1),XN(i,2),XN(i,3),num2str(i),'BackgroundColor','y') 
%     end
%     view(3)
    
    %% Newinds
%     I = [1,2,4,3,5,8,6,7];
    J = [1,3,7,9,19,21,25,27];
    INDS = zeros(27,1); %Indices
%     INDS(J) = iv(I); % Existing Corner nodes
    INDS(J) = iv(); % Existing Corner nodes
    
    NewNodes = (nmax+1:nmax+19)'; %Temporarly create
    
    %% Check if neighbors already have created same nodes
    % Temporarly add new coordinates to the list, new coordinates may
    % contain points that already exist
    X2 = [X;XN];
    % Those extra added points are found and stored as nodenumbers
    duplicateNodes = NewNodes( ismember(X2(NewNodes,:),X,'rows') );
    
    % The duplicates that we want to keep
    %     DupNodes = find(ismember(X,XN,'rows'));
    
    % First column stores indices to INDS second column stores node indices
    K = zeros(19,2);
    % We want to populate the noncorner nodes
    K(1:19,1) = [2,4,5,6,8,10,11,12,13,14,15,16,17,18,20,22,23,24,26]';
    % Loop over all 19 non corner nodes and look if the new coordinate
    % already exists in the domain. If it exists, add the index of the
    % existing node to the list, if it does not exist, increment the
    % nodenumber and add it to the liso.
    for i = 1:19
        ind = find(ismember(X,XN(i,:),'rows'),1);
        if ~isempty(ind)
            K(i,2)=ind;
        else
            K(i,2)=nmax+1;
            nmax=nmax+1;
        end
    end
    
    % Fill the local INDS list
    INDS(K(:,1)) = K(:,2);
    
    % The new global nodes in a connectivity matrix, ready to be appended to
    % the rest of the Connectivity matrix
    C2 = INDS(locnodes);
    
    % Delete the duplicate points from the new Coordinate matrix.
    X2(duplicateNodes,:) = [];
    % Replace coord matrix
    X = X2;
    
    % Add new element matrix to the rest
    nodes = [nodes;C2];
    
    dnodes = ismember(nodes,oldNodes(iel,:),'rows');
    nodes(dnodes,:)=[]; % Delete mother element, since it is being replaced
    NewElements(ic).old = iel;
    newElements = find(ismember(nodes,C2,'rows'));
    NewElements(ic).new = newElements;
    
    

    
    %% Hanging nodes
    EdgeNodes = find(ismember(X,XN(XNeind,:),'rows'));
    FaceNodes = find(ismember(X,XN(XNfind,:),'rows'));
    MidNode = find(ismember(X,XN(XNMind,:),'rows'));
    
    o.Element(iel).HangNodes.EdgeNodes = EdgeNodes;
    o.Element(iel).HangNodes.FaceNodes = FaceNodes;
    o.Element(iel).HangNodes.MidNode = MidNode;
    o.Element(iel).HangNodes.ParentElement = iel;
    
    
ic = ic+1;   
end



%% Viz
% xfigure(43434)
% fele = [6*iel-5;6*iel-4;6*iel-3;6*iel-2;6*iel-1;6*iel-0;];
% h.patch = patch(o.XC(o.Faces(fele(:),:)'),o.YC(o.Faces(fele(:),:)'),o.ZC(o.Faces(fele(:),:)'),'w','FaceColor','none');
% xlabel('X'); ylabel('Y'); zlabel('Z')
% axis equal tight
% 
% for i = 1:length(X)
%     text(X(i,1),X(i,2),X(i,3),num2str(i),'BackgroundColor','y')
% end
% view(3)


%% HangNodes handling
HangNodes = [o.Element.HangNodes];
ParentElement = [HangNodes.ParentElement];
EdgeNodes = [HangNodes.EdgeNodes];
FaceNodes = [HangNodes.FaceNodes];
MidNode = [HangNodes.MidNode];
HangNodesM = [ParentElement;EdgeNodes;FaceNodes;MidNode];
o.HangNodes = HangNodes;
% Hangnodes are ready to be returned
% Return HangNodes, HangNodesM

%% Update Face and edge matrices
Q = zeros(o.nele*6,4); %Stores faces
edges1 = zeros(o.nele*12,2); %stores edges

o.Element(o.nele).faces = []; %initialize o.Element.faces
o.Element(o.nele).edges = []; %initialize o.Element.edges

% I = [1,2,4,3,5,8,6,7];
% nodes = nodes(:,I);

o.nele = size(nodes,1);
lof = 1; loe = 1;
xc = X(:,1); yc = X(:,2); zc = X(:,3);
for iel = 1:o.nele
    upf = lof+5;
    upe = loe+11;
    %Faces
%     iface = [nodes(iel,[1,2,4,3]);...
%             nodes(iel,[5,6,8,7]);...
%             nodes(iel,[1,5,6,2]);...
%             nodes(iel,[2,6,8,4]);...
%             nodes(iel,[4,8,7,3]);...
%             nodes(iel,[3,7,5,1])];
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
%     I = [1,2,4,3,5,8,6,7]; %Numbering order (indices)
    iedges = [n(1),n(2);...
            n(2),n(4);...
            n(4),n(3);...
            n(3),n(1);...
            n(1),n(5);...
            n(3),n(7);...
            n(4),n(8);...
            n(2),n(6);...
            n(6),n(5);...
            n(5),n(7);...
            n(7),n(8);...
            n(8),n(6)];
    edges1(loe:upe,:) = iedges;
    o.Element(iel).edges = iedges;
    %counter
    lof = upf +1;
    loe = upe +1;
end

o.Connectivity = nodes;
o.Points = X;
o.XC = X(:,1);
o.YC = X(:,2);
o.ZC = X(:,3);

o.Faces = Q;
o.nnod = length(o.XC);
o.edges = edges1;

%% Viz mesh
% xfigure(43455)
% ele = 1:size(nodes,1);
% fele = [6*ele-5;6*ele-4;6*ele-3;6*ele-2;6*ele-1;6*ele-0;];
% h.patch = patch(o.XC(o.Faces(fele(:),:)'),o.YC(o.Faces(fele(:),:)'),o.ZC(o.Faces(fele(:),:)'),'w','FaceColor','none');
% xlabel('X'); ylabel('Y'); zlabel('Z')
% axis equal tight
% for i = 1:size(o.Points,1)
%     text(o.Points(i,1),o.Points(i,2),o.Points(i,3),num2str(i),'BackgroundColor','w')
% end
% view(3)
end
