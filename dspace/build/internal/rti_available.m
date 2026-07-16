function available = rti_available()
%RTI_AVAILABLE 探测 dSPACE RTI1202 (MicroLabBox) 块库是否可用
%
%   仅在 dSPACE 主机（Windows + RCP/HIL 软件）上为 true。
%   Linux 开发机上 build_rt_model 自动降级为 Inport/Outport 占位 IO。

available = false;
try
    if exist('rtilib1202', 'file') == 4 || ~isempty(which('rtilib1202'))
        load_system('rtilib1202');
        available = true;
    end
catch
    available = false;
end
end
