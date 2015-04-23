function h = vizMesh(o,varargin)
    %vizMesh Visualize the mesh in a figure
    %
    %   h = vizMesh()
    %   h = vizMesh(ele,figNumber,properties)
    %   ele is a list of elements to display
    %   h contains visualization properties
    %   figNumber sets the figure number
    %   
    %   Requirements: xfigure
    %
    %   properties:
    %   'NodeNumbers'       -   Display node numbers
    %   'ElementNumbers'    -   Display element numbers
    %
    %   Usage:
    %   o.vizMesh( ... )
    %   or
    %   vizMesh(T, ... )
    %   T is the Hex1Mesh class
    
    %
    ele = 1:size(o.Connectivity,1);
    if nargin > 1
        if isa(varargin{1},'double')
            ele = varargin{1};
        end
    end
    fn = NaN;
    if nargin > 2
        if isa(varargin{1},'double')
            fn = varargin{2};
        end
    end
    if isnan(fn)
        if exist('xfigure','file') == 2
            [h.fig,xf1] = xfigure;
        else
            RequiredFileMissing('xfigure', 'https://raw.githubuserconteno.com/cenmir/xfigure/master/xfigure.m')
            RequiredFileMissing('xfigure_KPF', 'https://raw.githubuserconteno.com/cenmir/xfigure/master/xfigure_KPF.m')
            [h.fig,xf1] = xfigure;
        end
    else
        if exist('xfigure','file') == 2
            [h.fig,xf1] = xfigure;
        else
            RequiredFileMissing('xfigure', 'https://raw.githubuserconteno.com/cenmir/xfigure/master/xfigure.m')
            RequiredFileMissing('xfigure_KPF', 'https://raw.githubuserconteno.com/cenmir/xfigure/master/xfigure_KPF.m')
            [h.fig,xf1] = xfigure;
        end
    end

    ele = ele(:);
    fele = [6*ele-5;6*ele-4;6*ele-3;6*ele-2;6*ele-1;6*ele-0;];

    h.patch = patch(o.XC(o.Faces(fele(:),:)'),o.YC(o.Faces(fele(:),:)'),o.ZC(o.Faces(fele(:),:)'),'w','FaceColor','none');
    xlabel('X'); ylabel('Y'); zlabel('Z')
    axis equal tight
    set(h.fig,'name','Hex1Mesh')
    view(-55,45)

    title(['Number of elements: ',num2str(o.nele),' Number of nodes:',num2str(o.nnod)])
    
    h = drawNodeNumbers(o,h,ele);
    
    if isenabled('NodeNumbers',varargin)
        if isfield(h,'NodeText')
            set(h.NodeText,'Visible','on')
        end
    end
    h = drawElementNumbers(o,h,ele);
    
    if isenabled('ElementNumbers',varargin)
        if isfield(h,'EleText')
            set(h.EleText,'Visible','on')
        end
    end
    
    
    helpText = xf1.uiTextHelp.String;
    helpText{end+1} = 'e - Toggle element numbers on or off';
    helpText{end+1} = 'n - Toggle node numbers on or off';        
    xf1.uiTextHelp.String = helpText;  
    xf1.uiTextHelp.Position = xf1.uiTextHelp.Position + [0,0,0,20];
    
    h.fig.KeyPressFcn = {@GKPFVizMesh,xf1,h};

end

function h = drawNodeNumbers(o,h,ele)
    if o.nele > 100
        warning('Cannot draw NodeNumbers, too many elements')
        return
    end
    h.NodeText = [];
    %                 unique(o.Connectivity(:),'stable')'
    %                 1:o.nnod
    nods = unique(o.Connectivity(ele,:))';
    for i = nods
        h.NodeText = [h.NodeText; text(o.XC(i),o.YC(i),o.ZC(i),num2str(i),'BackgroundColor','w') ];
    end
    set(h.NodeText,'Visible','off')
    
end

function h = drawElementNumbers(o,h,ele)
    if o.nele > 100
        warning('Cannot draw ElementNumbers, too many elements')
        return
    end

    h.EleText = [];
    for i = sort(ele)'
        xm = mean(o.XC(o.Connectivity(i,:)));
        ym = mean(o.YC(o.Connectivity(i,:)));
        zm = mean(o.ZC(o.Connectivity(i,:)));
        h.EleText = [h.EleText; text(xm,ym,zm,num2str(i),'BackgroundColor','y')];
    end
    set(h.EleText,'Visible','off')
end

function RequiredFileMissing(filename, RemoteDestination)
    %If we're going to download a whole bunch of files, it is better to set
    % RequiredFilesDir to be a global and not have to ask the user to
    % specify a destination folder for every file...
    GetABunchOfFiles = 1;
    
    if GetABunchOfFiles
        global RequiredFilesDir
    else
        RequiredFilesDir = [];
    end  
    
    warning([filename,' is missing!'])
    disp([filename,' is missing!'])
    disp(['Trying to download ',filename,' from ',RemoteDestination])
    
    
    if isempty(RequiredFilesDir)
        scriptPath = mfilename('class');
        [ScriptDir,~,~] = fileparts(scriptPath);
        DestDir = uigetdir(ScriptDir,['Select where to save ',filename,'. Make sure its either the script directory or a directory on the Path.']);
        if DestDir == 0
            error(['Failed to select folder, failed to install ',filename])
        end
        
        RequiredFilesDir = DestDir;
    end
    DestFile= [RequiredFilesDir,'/',filename];
    
    % Download the RemoteFile and save it to DestFile
    websave(DestFile,RemoteDestination);
    
    % Give up to 10 seconds for the file to show up, otherwise send error
    % message.
    tic
    while 1
        if exist('xfigure','file') == 2
            break
        end
        pause(0.1)
        t1 = toc;
        if t1 > 10
            error(['Failed to download ',filename,'! Timeouo.'])
        end
    end

    
end

function GKPFVizMesh(src,evnt,xfigure_This,h)
    %GKPF(src,evnt,xfigure_This,h)
    %   h is a struct that contains user defined graphical handles to be used
    %   in this function.
    %   xfigure_This is a handle to xfigure, returned by xfigure.

    %% Do not touch
    %Leave this be, we need the standard keypress function!
    xfigure_KPF(src, evnt, xfigure_This); %Do not touch
    
    
    %% Add own keys below

    if strcmpi(evnt.Character,'e')
        if any(strcmpi(get(h.EleText,'Visible'),'off'))
            set(h.EleText,'Visible','on')
        else
            set(h.EleText,'Visible','off')
        end
    end

    if strcmpi(evnt.Character,'n')
        if any(strcmpi(get(h.NodeText,'Visible'),'off'))
            set(h.NodeText,'Visible','on')
        else
            set(h.NodeText,'Visible','off')
        end
    end


end



