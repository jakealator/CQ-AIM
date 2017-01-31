% Gread.m
% Created by Peter Monk
% Mesh=Gread(file)
% Prototype to read a Gmsh file into a matlab data structure called Mesh
%
% It is important that "physical" tags have been used to generate the mesh
% otherwise boundary region information will be lost
%
% Mesh data structure (for a 2D mesh)
%
%         NPL: Number of physical labels
%         Pnames: {3x3 cell} labels (number, index, text)
%         NoN: number of nodes
%         nodes: NoN x 3 array of nodes in 3D
%         Ne: number of edges
%         Nf: number of triangles
%         Nt: number of tets
%         Nv: number of vertices (labeled)
%         edge: Ne x 2 array of end points of each edge
%         edgetag: Ne x ? array of tags (the first is physical)
%         faces: Nf x 3 array of vertices in faces
%         facetag: Nf x ? array of tags for triangles
%

function Mesh=Gread(file)

meshfile=[file,'.msh'];
NPL=0;

[fid message] = fopen(meshfile,'r');
if fid==-1
    disp('Unable to open Gmesh file :-(')
    disp('Heres what the computer said:')
    disp(message)
    error('Fatal error, :-(')
else
    while feof(fid)==0
        L=fgetl(fid);
        %meshformat
        if strcmp(L,'$MeshFormat')
            inp=fscanf(fid,'%e %i %i');
            disp(sprintf('Gmsh Version = %3.1f',inp(1)));
            L=fgetl(fid);
        end
        % pysical names
        if strcmp(L,'$PhysicalNames')
            Mesh.NPL= fscanf(fid,'%i',1);
            if Mesh.NPL>0, NPL=1;,end
            Mesh.Pnames=cell(Mesh.NPL,3);
            for j=1:Mesh.NPL
                Mesh.Pnames{j,1}=fscanf(fid,'%i',1);
                Mesh.Pnames{j,2}=fscanf(fid,'%i',1);
                tmp=textscan(fid,'%s%*[^\n]',1);
                label=tmp{1}{1};
                label=label(2:length(label)-1);
                if Mesh.Pnames{j,1}==1
                    disp([sprintf('Edges with Gmsh index %i are %a',...
                        Mesh.Pnames{j,2}),label])
                elseif Mesh.Pnames{j,1}==2
                    disp([sprintf('Triangles with Gmsh index %i are %a',...
                        Mesh.Pnames{j,2}),label])
                elseif Mesh.Pnames{j,1}==3
                    disp([sprintf('Regions with Gmsh index %i are %a',...
                        Mesh.Pnames{j,2}),label])
                else
                    disp('Unknown labeled entity')
                end
                Mesh.Pnames{j,3}=label;
            end
            L=fgetl(fid); % get the empty line?
            L=fgetl(fid);
        end
        %nodes
        if strcmp(L,'$Nodes')
            Mesh.NoN= fscanf(fid,'%i',1);
            tmp=fscanf(fid,'%g',4*Mesh.NoN);
            tmp=reshape(tmp,4,Mesh.NoN);
            tmp=tmp';
            Mesh.nodes=tmp(:,2:4);
            L=fgetl(fid);
            L=fgetl(fid);
        end
        % elements
        if strcmp(L,'$Elements')
            NE= fscanf(fid,'%i',1);
            disp(sprintf('Reading %i elements',NE))
            Mesh.Ne=0;
            Mesh.Nf=0;
            Mesh.Nt=0;
            Mesh.Nv=0;
            for j=1:NE
                idno=fscanf(fid,'%i',1);
                elt=fscanf(fid,'%i',1);
                not=fscanf(fid,'%i',1);
                tags=fscanf(fid,'%i',not);
                if elt==1
                    % edge
                    verts=fscanf(fid,'%i',2);
                    Mesh.Ne=Mesh.Ne+1;
                    Mesh.edge(Mesh.Ne,:)=verts'; 
                    Mesh.edgetag(Mesh.Ne,:)=tags;
                 elseif elt==4
                    % tet
                    verts=fscanf(fid,'%i',4);
                    Mesh.Nt=Mesh.Nt+1;
                    Mesh.tets(Mesh.Nt,:)=verts';
                    Mesh.tetstag(Mesh.Nt,:)=tags;
                elseif elt==2
                    % triangle
                    verts=fscanf(fid,'%i',3);
                    Mesh.Nf=Mesh.Nf+1;
                    Mesh.faces(Mesh.Nf,:)=verts';
                    Mesh.facetag(Mesh.Nf,:)=tags;
               elseif elt==15
                    % nodes nodes
                    nodes=fscanf(fid,'%i',1);
                    Mesh.Nv=Mesh.Nv+1;
                    Mesh.verts(Mesh.Nv,:)=nodes;
                    Mesh.vertstag(Mesh.Nv,:)=tags;
                else
                    disp(sprintf('Element type = %i',elt))
                    error('Error: Element type not implemented')
                end
            end
            L=fgetl(fid);
            L=fgetl(fid);
        end
    end
end
if NPL==0
    disp('Warning: no physical labels in this data')
end
