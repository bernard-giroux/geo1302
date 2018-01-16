classdef MaillageTriangulaire < Maillage
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        noeuds
        triangles
    end
    properties (SetAccess=private)
        nom
    end

    methods
        function obj = MaillageTriangulaire(nom, n, t)
            if nargin > 0
                obj.nom = nom;
                obj.noeuds = n;
                obj.triangles = t;
            end
        end
        function obj = set.noeuds(obj, n)
            if isnumeric(n)
                if size(n,2)~=3
                    error('n doit etre de dimension nNoeuds x 3')
                end
                obj.noeuds = n;
            else
                error('Les arguments doivent etre numeriques')
            end
        end
        function obj = set.triangles(obj, t)
            if isnumeric(t)
                if size(t,2)~=3
                    error('n doit etre de dimension nNoeuds x 3')
                end
                obj.triangles = t;
            else
                error('Les arguments doivent etre numeriques')
            end
        end
        function obj = set.nom(obj, n)
            if ischar(n)
                obj.nom = n;
            else
                error('Le nom doit etre une chaine de caracteres')
            end
        end
        function obj = setNom(obj, n)
            obj.nom = n;
        end

        function plot(obj, varargin)
            patch('faces',obj.triangles,'vertices',obj.noeuds,varargin{:});
            view(3)
            axis vis3d
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]')
            camlight
        end

        function lim = getLimites(obj)
            lim = [min(obj.noeuds) max(obj.noeuds)];
        end
        function n = getNombreDeNoeuds(obj)
            n = size(obj.noeuds,1);
        end
        function n = getNombreDeCellules(obj)
            n = size(obj.triangles,1);
        end
    end

end
