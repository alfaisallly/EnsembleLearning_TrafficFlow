classdef Config

    methods(Static)
        
        function [x] = get_names()
            x = {'smalltest','210E'};
        end
        
        function [x] = get(configname)
            x = struct( 'xls_file' , '' , ...
                        'xml_file' , '' , ...
                        'sim_dt' , nan , ...
                        'vds_use_template' , [] , ...
                        'good_days' , [] , ...
                        'model_day' , nan );
                    
            switch configname
                
                case 'smalltest'
                    x.xml_file  = '_smalltest.xml';
                    x.xls_file  = '';
                    x.sim_dt    = 5;
                    x.vds_use_template  = [];
                    x.model_day = nan;
                    x.good_days = [];
                    
                case '210E'
                    x.xml_file  = '210E_generated_x43_sensors_norepeats.xml';
                    x.xls_file  = 'I210EB_Data_x43.xlsx';
                    x.sim_dt    = 4;
                    
                    i=1;
                    x.vds_use_template(i).link_id  = 135273082;
                    x.vds_use_template(i).id       = 716563;
                    x.vds_use_template(i).type     = 'OR';
                    x.vds_use_template(i).location = 'Colorado Blvd';
			 			
                    i=i+1;
                    x.vds_use_template(i).link_id  = 756090723;
                    x.vds_use_template(i).id       = 770172;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = '210';

                    i=i+1;
                    x.vds_use_template(i).link_id  = 126499392;
                    x.vds_use_template(i).id       = 8800004;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = '710 (Corson) 34.1499 -118.1549';

                    i=i+1;
                    x.vds_use_template(i).link_id  = -873337960;
                    x.vds_use_template(i).id       = 9900001;
                    x.vds_use_template(i).type     = 'OR';
                    x.vds_use_template(i).location = 'Fair Oaks HOV entrance';

                    i=i+1;
                    x.vds_use_template(i).link_id  = 128793003;
                    x.vds_use_template(i).id       = 716589;
                    x.vds_use_template(i).type     = 'OR';
                    x.vds_use_template(i).location = 'N Hill Ave';

                    i=i+1;
                    x.vds_use_template(i).link_id  = -126499580;
                    x.vds_use_template(i).id       = 763908;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = 'N San Gabriel Blvd	';
			 		
                    i=i+1;
                    x.vds_use_template(i).link_id  = -126537863;
                    x.vds_use_template(i).id       = 761167;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = 'S Mountain Ave	';
                    
                    i=i+1;
                    x.vds_use_template(i).link_id  = -24026218;
                    x.vds_use_template(i).id       = 8800003;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = 'Buena Vista St 34.1354 -117.9863';

                    i=i+1;
                    x.vds_use_template(i).link_id  = -781754904;
                    x.vds_use_template(i).id       = 769705;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = '605 SB';

                    i=i+1;
                    x.vds_use_template(i).link_id  = -24558231;
                    x.vds_use_template(i).id       = 8800001;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = '605 NB (Mnt Olive) 34.1351 -117.9612';

                    i=i+1;
                    x.vds_use_template(i).link_id  = 24162804;
                    x.vds_use_template(i).id       = 769774;
                    x.vds_use_template(i).type     = 'OR';
                    x.vds_use_template(i).location = '605 NB';

                    i=i+1;
                    x.vds_use_template(i).link_id  = -781749995;
                    x.vds_use_template(i).id       = 774990;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = 'N Irwindale Ave';

                    i=i+1;
                    x.vds_use_template(i).link_id  = -24029866;
                    x.vds_use_template(i).id       = 8800002;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = 'N Vernon Ave 34.1282 -117.9201';

                    i=i+1;
                    x.vds_use_template(i).link_id  = -24558314;
                    x.vds_use_template(i).id       = 718213;
                    x.vds_use_template(i).type     = 'OR';
                    x.vds_use_template(i).location = 'MOUNT OLIVE DR / 605 34.13449	-117.958916';

                    %**************************************************
                    %The following links/vdss are selected from
                    %I210EB_Data_x43 
                    i=i+1;
                    x.vds_use_template(i).link_id  = -121240605;
                    x.vds_use_template(i).id       = 717594;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = 'FR N. Figueroa St';
                    
                    i=i+1;
                    x.vds_use_template(i).link_id  = -128793036;
                    x.vds_use_template(i).id       = 717671;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = 'FR N Santa Anita Ave';
                    
                    i=i+1;
                    x.vds_use_template(i).link_id  = -24026240;
                    x.vds_use_template(i).id       = 761154;
                    x.vds_use_template(i).type     = 'FR';
                    x.vds_use_template(i).location = 'FR S Myrtle Ave';
                    
                    i=i+1;
                    x.vds_use_template(i).link_id  = 24031575;
                    x.vds_use_template(i).id       = 717679;
                    x.vds_use_template(i).type     = 'OR';
                    x.vds_use_template(i).location = 'OR N Azusa Ave';
                    
                    %**************************************************
                    
                    x.model_day = datenum(2014,10,13);
                    
                    % REMOVED THIS DAY FROM x.good_days
                    % datenum('8-Nov-14'),...  
                    
                    x.good_days = [ datenum('1-Oct-14'),...
                                    datenum('2-Oct-14'),...
                                    datenum('3-Oct-14'),...
                                    datenum('11-Oct-14'),...
                                    datenum('12-Oct-14'),...
                                    datenum('13-Oct-14'),...
                                    datenum('14-Oct-14'),...
                                    datenum('15-Oct-14'),...
                                    datenum('17-Oct-14'),...
                                    datenum('18-Oct-14'),...
                                    datenum('19-Oct-14'),...
                                    datenum('20-Oct-14'),...
                                    datenum('21-Oct-14'),...
                                    datenum('22-Oct-14'),...
                                    datenum('23-Oct-14'),...
                                    datenum('29-Oct-14'),...
                                    datenum('30-Oct-14'),...
                                    datenum('31-Oct-14'),...
                                    datenum('1-Nov-14'),...
                                    datenum('2-Nov-14'),...
                                    datenum('7-Nov-14'),...
                                    datenum('12-Nov-14'),...
                                    datenum('13-Nov-14'),...
                                    datenum('14-Nov-14'),...
                                    datenum('15-Nov-14'),...
                                    datenum('17-Nov-14'),...
                                    datenum('18-Nov-14'),...
                                    datenum('25-Nov-14'),...
                                    datenum('15-Dec-14'),...
                                    datenum('16-Dec-14'),...
                                    datenum('17-Dec-14'),...
                                    datenum('18-Dec-14'),...
                                    datenum('19-Dec-14'),...
                                    datenum('20-Dec-14'),...
                                    datenum('31-Dec-14')]; 
                                
%                     x.good_days = x.good_days(6:7);
                                
            end
            
            if isempty(x.vds_use_template)
                x.vds_use_template_ids = [];
            else
                x.vds_use_template_ids = [x.vds_use_template.id];
            end
            
            
        end
        
    end
    
end
