classdef Folder
    
    methods(Static)
                
        function [x] = root()
            here = fileparts(mfilename('fullpath'));
            x = fileparts(here);
        end
        
        function [x] = source()
            x = fullfile(Folder.root,'src');
        end
        
        function [x] = reports()
            x = fullfile(Folder.root,'reports');
        end
        
        function [x] = data()
            x = fullfile(Folder.root,'data');
        end
        
        function [x] = config()
            x = fullfile(Folder.root,'config');
        end
        
        function [x] = tests()
            x = fullfile(Folder.source,'tests');
        end
        
        function [x] = temp()
            x = fullfile(Folder.root,'temp');
        end
        
        function [x] = objects()
            x = fullfile(Folder.root,'obj');
        end
        
    end
    
end

