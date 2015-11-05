function this = writeCcode_sensi(this,svar,fid)
% writeCcode_sensi is a wrapper for writeCcode which loops over parameters and reduces
% overhead by check nonzero values
%
% Parameters:
%  svar: symbolic variable which is later on passed to ggode @type symbolic
%  fid: file id in which the final expression is written @type fileid
% 
% Return values:
%  this: model definition object @type amimodel
      
np = this.np;

if(strcmp(svar,'drvaldp'))
    nonzero = this.sym.drootpdp~=0;
    if(any(any(nonzero)))
        for ip=1:np
            if(any(nonzero(:,ip)))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                this.writeCcode('drvaldp', fid, ip);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(svar,'s2root'))
    nonzero = this.sym.s2root~=0;
    if(any(any(any(nonzero))))
        for ip=1:np
            if(any(any(nonzero(:,:,ip))))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                fprintf(fid,'    sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                fprintf(fid,['      switch(jp) {\n']);
                for jp=1:np
                    if(any(nonzero(:,jp,ip)))
                        fprintf(fid,['          case ' num2str(jp-1) ': {\n']);
                        this.writeCcode('s2root', fid, ip, jp);
                        fprintf(fid,'\n');
                        fprintf(fid,'          } break;\n\n');
                    end
                end
                fprintf(fid,'\n      }\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(svar,'s2rootval'))
    nonzero = this.sym.s2rootval~=0;
    if(any(any(any(nonzero))))
        for ip=1:np
            if(any(any(nonzero(:,:,ip))))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                fprintf(fid,'  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                fprintf(fid,['      switch(jp) {\n']);
                for jp=1:np
                    if(any(nonzero(:,jp,ip)))
                        fprintf(fid,['          case ' num2str(jp-1) ': {\n']);
                        this.writeCcode('s2rootval', fid, ip, jp);
                        fprintf(fid,'\n');
                        fprintf(fid,'           } break;\n\n');
                    end
                end
                fprintf(fid,'\n      }\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(svar,'ideltadisc'))
    nonzero = this.sym.(svar) ~=0;
    if(any(nonzero))
        for ip=1:np
            if(nonzero(ip))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                this.writeCcode(svar, fid, ip);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
else 
    nonzero = this.sym.(svar) ~=0;
    if(any(any(nonzero)))
        for ip=1:np
            if(any(nonzero(:,ip)))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                if(strcmp(svar,'sroot') || strcmp(svar,'srootval') || strcmp(svar,'sdeltadisc'))
                    fprintf(fid,'  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                end
                this.writeCcode(svar, fid, ip);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
end
end