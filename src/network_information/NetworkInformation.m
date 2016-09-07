classdef NetworkInformation < handle

    properties
        configfile
        scenario    @ScenarioPtr
        link_ids 
        node_ids
        link_lengths_meters
        vds2linkid
    end
    
    methods(Access=public)
        
        function [this] = NetworkInformation(configfile)
            this.configfile = configfile;
            this.scenario = ScenarioPtr;
            this.scenario.load(configfile);
            this.link_ids = Utils.row_vector( this.scenario.get_link_ids );
            this.node_ids = Utils.row_vector(this.scenario.get_node_ids);
            this.link_lengths_meters = Utils.row_vector(this.scenario.get_link_lengths('si'));
            T = this.scenario.get_sensor_table;
            
            if ~isempty(T)
                is_attached = logical(T.is_attached);               
                this.vds2linkid = [T.vds(is_attached) T.link_id(is_attached)];
                this.vds2linkid = this.vds2linkid';
            end
        end
        
        function [x]=get_source_link_ids(this)
            x = this.link_ids(this.scenario.is_source_link);
        end
        
        function [x] = get_internal_link_ids(this)
            x = this.scenario.get_link_ids(~this.scenario.is_source_link & ~this.scenario.is_sink_link);
            x = Utils.row_vector(x);
        end
        
        function [x]=get_node_io_for_multi_output_nodes(this)
            node_io = this.scenario.get_node_io(this.node_ids);
            multi_out = false(1,length(node_io));
            for i=1:length(node_io)
                multi_out(i) = length(node_io(i).link_out) > 1;
            end
            x = node_io(multi_out);
        end
        
        function [x]=get_all_vds(this)
            if ~isempty(this.vds2linkid)
                x = this.vds2linkid(1,:);
            else
                x = [];
            end
        end
        
        function [lids] = get_link_ids_for_vds(this,vds)
            if ~isempty(this.vds2linkid)
                lids = this.vds2linkid(2,Utils.index_into(vds,this.vds2linkid(1,:)));
            else
                lids = [];
            end
        end
        
        function [vdss] = get_vds_for_link_ids(this,lids)
            if ~isempty(this.vds2linkid)
                vdss = nan(1,length(lids));
                [isthere,ind] = ismember(lids,this.vds2linkid(2,:));
                vdss(isthere) = this.vds2linkid(1,ind(isthere));
            else
                vdss = [];
            end
        end
        
        function [L] = get_lengths_for_linkids_in_km(this,linkids)
            link_lengths = this.scenario.get_link_lengths('si')/1000;
            L = link_lengths(Utils.index_into(linkids,this.link_ids));
            L = Utils.row_vector(L);
        end
        
        % get next downstream (direction>0) links or upstream (direction<0) links
        function [x] = get_next_links(this,thislink,direction)
            A = this.scenario.link_id_begin_end;
            ind = A(:,1)==thislink;
            if ~any(ind)
                error('link id not found')
            end
            
            if direction>0      % downstream links
                x = A(A(:,2)==A(ind,3),1);
            else                % upstream links
                x = A(A(:,3)==A(ind,2),1);
            end 
        end
        
        function [b,e] = get_begin_end_node_for_link(this,thislink)
            ind = this.scenario.link_id_begin_end(:,1)==thislink;
            if ~any(ind)
                error('~any(ind)')
            end
            b = this.scenario.link_id_begin_end(ind,2);
            e = this.scenario.link_id_begin_end(ind,3);
        end
           
        function [in,out] = get_in_out_links_for_node(this,this_node)
            
            if ~ismember(this_node,this.node_ids)
                error('~ismember(this_node,this.node_ids)')
            end
            
            ind = this.scenario.link_id_begin_end(:,3)==this_node;
            in = this.scenario.link_id_begin_end(ind,1)';
            
            ind = this.scenario.link_id_begin_end(:,2)==this_node;
            out = this.scenario.link_id_begin_end(ind,1)';

        end
        
        function [x] = get_link_is_terminal(this,link_id)
            ind = this.scenario.link_id_begin_end(:,1)==link_id;
            if ~any(ind)
                error('~any(ind)')
            end
            begin_end_nodes = this.scenario.link_id_begin_end(ind,[2 3]);
            node_io = this.scenario.get_node_io(begin_end_nodes);
            x = isempty(node_io(1).link_in) || isempty(node_io(2).link_out);
        end
        
        
    end
    
end
