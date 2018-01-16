classdef (Abstract) Maillage
    %MAILLAGE Classe abstraite pour representer des maillages
   methods
     lim = getLimites(obj)
     nn = getNombreDeNoeuds(obj)
     nc = getNombreDeCellules(obj)
   end
end

