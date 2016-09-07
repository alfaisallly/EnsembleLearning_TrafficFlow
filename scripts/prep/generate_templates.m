function []=generate_templates()
% Create pems-like structures and mat files for all onramps and offramps.
% Temlate flows are taken from the configuration file.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(fullfile(root,'src')))

name = '210E';
config = Config.get(name);
ni = ObjectFactory.network_information(config.xml_file);

[xls_ors,xls_frs] = get_xls_links(fullfile(Folder.config,config.xls_file));
[or,fr]=read_xls(fullfile(Folder.config,config.xls_file));
source_sink_link_ids = setdiff(ni.link_ids,ni.get_internal_link_ids);
   
template = [];

% generate data for boundaries from beats run
beats = BeatsSimulation;
beats.load_scenario(fullfile(Folder.config,config.xml_file));
beats.run_beats(struct('SIM_DT',config.sim_dt));
fwy = beats.scenario_ptr.get_freeway_structure;

%  put upstream it into template
up_ml_link = fwy.linear_fwy_link_ids(1);
ind = fwy.linear_fwy_sensors.link_ids==up_ml_link;
vdss = fwy.linear_fwy_sensors.sensor_vds(ind);
for j=1:length(vdss)
    template = [template ...
        struct('vds',vdss(j), ...
        'link_id',up_ml_link, ...
        'time_sec',0:300:86100, ...
        'flow_vph',beats.get_output_for_link_id(up_ml_link).flw_out_vph' ) ];
end

%  put downstream it into template
dn_ml_link = fwy.linear_fwy_link_ids(end);
ind = fwy.linear_fwy_sensors.link_ids==dn_ml_link;
vdss = fwy.linear_fwy_sensors.sensor_vds(ind);
for j=1:length(vdss)
    template = [template ...
        struct('vds',vdss(j), ...
        'link_id',dn_ml_link, ...
        'time_sec',0:300:86100, ...
        'flow_vph',beats.get_output_for_link_id(dn_ml_link).flw_out_vph' ) ];
end

clear beats

for i = 1:length(source_sink_link_ids)
        
    link_id = source_sink_link_ids(i);
    
    is_or = ismember(link_id,xls_ors);
    is_fr = ismember(link_id,xls_frs);
    if ~xor(is_or,is_fr)
        warning('source or sink with no template. link_id = %d',link_id)
        continue
    end
    
    flow_vph = [];
    if is_or
        flow_vph = or(xls_ors==link_id,:);
    end
    if is_fr
        flow_vph = fr(xls_frs==link_id,:);
    end
      
    ind = ismember(ni.vds2linkid(2,:),link_id);
    if ~any(ind)
        error('~any(ind)')
    end
    
    vdss = ni.vds2linkid(1,ind);
  
    for j=1:length(vdss)
        template = [template ...
            struct('vds',vdss(j), ...
            'link_id',link_id, ...
            'time_sec',0:300:86100, ...
            'flow_vph',flow_vph ) ];
    end
                        
end

save(fullfile(Folder.data,sprintf('%s_template',name)),'template')

function [or,fr]=read_xls(xls_file)

or_nom = read_data(xls_file,'On-Ramp_CollectedFlows');
or_knb = read_data(xls_file,'On-Ramp_Knobs');
fr_nom = read_data(xls_file,'Off-Ramp_CollectedFlows');
fr_knb = read_data(xls_file,'Off-Ramp_Knobs');

or = or_nom.*or_knb;
fr = fr_nom.*fr_knb;

function [or_id,fr_id] = get_xls_links(xls_file)
[~,~,X] = xlsread(xls_file,'Configuration');
fr_id = cellfun(@(x)x,X(2:end,strcmp(X(1,:),'Off-Ramp ID')));
or_id = cellfun(@(x)x,X(2:end,strcmp(X(1,:),'On-Ramp ID')));

function [data,X,start_ind,end_ind]=read_data(xls_file,sheetname)
[~,~,X] = xlsread(xls_file,sheetname);
start_ind = find(cellfun(@(x)ischar(x),X(1,:)),1,'last')+1;
end_ind = length(X(1,:));
m = end_ind-start_ind+1;
n = size(X,1)-1;
x = X(2:end,start_ind:end_ind);
data = reshape([x{:}],n,m);

