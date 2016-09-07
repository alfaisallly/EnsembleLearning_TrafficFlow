classdef DataProfile

    properties
        id 
        day   
		time
        flw_in_vph 
		flw_out_vph		
		dty_vpk 
    end
	
    methods( Access = public )

        function [this] = DataProfile(id,day,time,flw_in_vph,flw_out_vph,dty_vpk)
        
            if nargin==0
                this.id = nan;
                this.day = nan;
                this.time = [];
                this.flw_in_vph = [];
                this.flw_out_vph = [];
                this.dty_vpk = [];
                return;
            end
            
            if isempty(time) 
                error('empty data')
            end
            
            n = length(time);
                        
            if ~isempty(flw_out_vph) && n-1~=length(flw_out_vph)
                error('bad dimensions')
            end   
            
            if ~isempty(dty_vpk) && n~=length(dty_vpk)
                error('bad dimensions')
            end    
            
            if ~isempty(flw_in_vph) && n-1~=length(flw_in_vph)
                error('bad dimensions')
            end
            
            this.id = id;
            this.day = day;
			this.time = Utils.row_vector(time);
			this.flw_in_vph = Utils.row_vector(flw_in_vph);
			this.flw_out_vph = Utils.row_vector(flw_out_vph);
			this.dty_vpk = Utils.row_vector(dty_vpk);
		end
		
	end

end

