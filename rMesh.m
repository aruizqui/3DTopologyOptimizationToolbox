classdef rMesh < handle
    % Generate or Import TET4 Mesh and store as node coordinates array and
    % element connectivity array. Can set elements to be removed (unless a
    % load is applied on them) and elements that cannot be removed in the
    % optimization. 
    % v1 Feb 2021 - V. Cholvi: Inital version Toolbox
    % v2 Dec 2024 - A. Ruiz: Script adapted to import either TET4 & HEX8
    %               New/Modified functions:
    %               --> importData2(ms.fileName)
    %               --> plot(obj, varargin)
    %               --> hexasurf(obj, varargin)
    
    properties
        elemConn % Element Connectivity Array: Lists Nodes in Each Element 
        nodeCoord % Node Coordinates: Lists X, Y, Z coords. for each node
        zeroElements    % Elements that have to be removed 
        mS              % Mesh Settings
        nonDesign       % List of elements that cannot be removed 
    end
    
    methods 
        function obj = rMesh(mS)
            obj.mS = mS;
            if isfield(mS, 'import') && mS.import
                [obj.elemConn, obj.nodeCoord] = importData2(mS.fileName);
                %obj.elemConn = importdata("conectivity.txt");
                %obj.nodeCoord = importdata("coordinates.txt");
            else
            [obj.elemConn, obj.nodeCoord] = tet4RelaxationMesh(mS);
            end
            obj.zeroElements = [];
            obj.nonDesign = [];
        end
        
        function nE = numElem(obj)
            % Returns Number of Elements
            nE = size(obj.elemConn, 1);
        end
        
        function toE = typeofElement(obj)
            %Return Number of Nodes per Element
            %CTETRA ==4; HEX8==8 
            toE = size(obj.elemConn,2);
        end

        function nN = numNodes(obj)                 
            % Returns Number of Nodes
            nN = size(obj.nodeCoord, 1);
        end
        
        function xCoord = X(obj, varargin)                    
            % Returns Nodal X-Coordinates
            x = obj.nodeCoord(:,1);
            if isempty(varargin) || isempty(varargin{1})
                xCoord = x;
            else 
                xCoord = x(varargin{1});
            end
        end
        
        function yCoord = Y(obj, varargin)                    
            % Returns Nodal Y-Coordinates
            y = obj.nodeCoord(:, 2);
            if isempty(varargin) || isempty(varargin{1})
                yCoord = y;
            else 
                yCoord = y(varargin{1});
            end
        end
        
        function zCoord = Z(obj, varargin)                    
            % Returns Nodal Z-Coordinates
            z = obj.nodeCoord(:,3);
            if isempty(varargin) || isempty(varargin{1})
                zCoord = z;
            else 
                zCoord = z(varargin{1});
            end
        end
        
        % Remove elements that contain these nodes
        function removeElements(obj, condition) 
            affectedNodes = find(condition);
            elems = false(obj.numElem,1);
            
            for i = 1:size(obj.elemConn, 2)
                elems = elems | ismember(obj.elemConn(:,i), affectedNodes);
            end
            
            obj.zeroElements = unique([obj.zeroElements; find(elems)]);
        end
        
        function nonDesignElements(obj, condition)
            % Elements that cannot be removed in the optimization
            affectedNodes = find(condition);
            elems = false(obj.numElem, 1);

            for i = 1:size(obj.elemConn, 2) 
                elems = elems | ismember(obj.elemConn(:,i), affectedNodes);
            end
            elems = sum(obj.elemConn(:,1) == affectedNodes', 2) | ...
                    sum(obj.elemConn(:,2) == affectedNodes', 2) | ...
                    sum(obj.elemConn(:,3) == affectedNodes', 2) | ...
                    sum(obj.elemConn(:,4) == affectedNodes', 2);
                
            obj.nonDesign = unique([obj.nonDesign; find(elems)]);
        end
        
        function plot(obj, varargin)                
        % Mesh Plotting. First Argument(optional) is array of values
        % to color differently for each element. Second
        % Argument(optional) is logical indexing of elements to be 
        % plotted
            if obj.typeofElement == 4
                %CTETRA
                disp('Mesh built with CTETRA elements')
                if length(varargin) > 1
                    elemsToPlot = setdiff(find(varargin{2}), obj.zeroElements);
                else
                    elemsToPlot = setdiff(1:obj.numElem, obj.zeroElements);                
                end
            
                if ~isempty(varargin) > 0 && (~isempty(varargin{1}))
                    colorValue = varargin{1};
                    obj.tetrasurf(elemsToPlot, colorValue)
                else
                    obj.tetrasurf(elemsToPlot, [0.8   0.8   1])
                end
                daspect([1 1 1])
            %HEX8
            elseif obj.typeofElement == 8
                disp('Mesh built with HEX8 elements')

                if length(varargin) > 1
                    elemsToPlot = setdiff(find(varargin{2}), obj.zeroElements);
                else
                    elemsToPlot = setdiff(1:obj.numElem, obj.zeroElements);
                end
                
                if ~isempty(varargin) && ~isempty(varargin{1})
                    colorValue = varargin{1};
                    obj.hexsurf(elemsToPlot, colorValue);
                else
                    obj.hexsurf(elemsToPlot, [0.8, 0.8, 1]);
                end
                daspect([1 1 1])
                
            else
                error('Unsupported element type');
            end
            
        end
        
        function tetrasurf(obj, etp, C)
            % Use Matlab Trisurf triangle plotting to plot tetrahedra
            % obj.tetrasurf(etp, C) etp: elements to plot, C: color/value
             if length(C) <= 6 && ~ (all(C == 1) && length(C) ==1)
                trisurf(obj.elemConn(etp,1:3),obj.X,obj.Y,obj.Z, ...
                    'Facecolor',C)
                hold on
                trisurf(obj.elemConn(etp,2:end), obj.X, obj.Y, obj.Z, ...
                    'Facecolor',C)
                trisurf(obj.elemConn(etp,[1 3 4]), obj.X, obj.Y, obj.Z, ...
                    'Facecolor',C)
                trisurf(obj.elemConn(etp,[1 2 4]), obj.X, obj.Y, obj.Z, ...
                    'Facecolor',C)
            else 

                trisurf(obj.elemConn(etp,1:3),obj.X,obj.Y,obj.Z,C(etp))
                hold on
                trisurf(obj.elemConn(etp,2:end),obj.X,obj.Y,obj.Z,C(etp))
                trisurf(obj.elemConn(etp,[1 3 4]),obj.X,obj.Y,obj.Z,C(etp))
                trisurf(obj.elemConn(etp,[1 2 4]),obj.X,obj.Y,obj.Z,C(etp))
            end
            axis tight
            
        end
        
        function hexsurf(obj, etp, C)
        % Use of `patch` to plot HEX8 elements.
        % etp: index of elements to plot.
        % C: color for each element.
            faces = [
                1, 2, 3, 4;
                5, 6, 7, 8;
                1, 4, 8, 5;
                2, 3, 7, 6;
                1, 2, 6, 5;
                4, 3, 7, 8;
                ];
            
            clf; %Clear previous plot
            axis tight
            view(3);
            hold on;
            
            if length(C)==1||length(C) ==3
               C_plot = repmat(C,length(etp),1); 
            else
                C_plot = C(etp);
            end

            for idx = 1:length(etp)
                elemidx = etp(idx);% Element index
                % Nodes of the element and coordinates
                elemNodes = obj.elemConn(elemidx, :);  
                elemCoords = obj.nodeCoord(elemNodes,:); 

                for j = 1:size(faces,1)
                    % Coord nodes of the face j
                    faceNodes = elemCoords(faces(j,:),:); 
                    if length(C) >1% obj.numElem
                    patch('Faces', [1 2 3 4],'Vertices',faceNodes,...
                        'FaceColor','flat', ... 
                        'FaceVertexCData', C_plot(idx),...
                        'EdgeColor', 'black',...%Edge color
                        'FaceAlpha',0.8) %Transparency
                    else
                    patch('Faces', [1 2 3 4],'Vertices',faceNodes,...
                        'FaceColor',C_plot(idx,:), ... 
                        'EdgeColor', 'black',...%Edge color
                        'FaceAlpha',0.8) %Transparency
                    end
                end
            end
            
            hold off;
        end
        
        function plotNonDesign(obj)
           if obj.typeofElement ==4
               obj.tetrasurf(obj.nonDesign, 'cyan')
           elseif obj.typeofElement ==8
               obj.hexsurf(obj.nonDesign, 'cyan')
           end
        end
    end
end