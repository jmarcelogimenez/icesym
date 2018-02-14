function [theta, time, hc, hr, dQhc, dQhr, xb, xbdot, ...
	  m_in, m_out, vol, vdot, m_fuel, m_air, m_res, ...
	  dQ_ht, dQ_chem] = read_cylextra(filename)
%
%
%

[bas, output] = system(["wc " filename]);
nlines = str2num(output(1:9));

if rem(nlines,4) != 0
  disp([" ERROR IN FILE " filename])
  error(" INCORRECT NUMBER OF LINES ")
end
nblocks = nlines/4;

old_cycle = 0;

fid = fopen(filename);
for iblock=1:nblocks
  # Reads line 1
  cycle = fscanf(fid, "%d", 1);
  if cycle != old_cycle
     ii = 0;
     old_cycle = cycle;
  end
  ii += 1;
  theta(ii,cycle) = fscanf(fid, "%f", 1);
  time(ii,cycle)  = fscanf(fid, "%f", 1);

  # Reads line 2
  hc(ii,cycle) = fscanf(fid, "%f", 1);
  hr(ii,cycle) = fscanf(fid, "%f", 1);
  dQhc(ii,cycle) = fscanf(fid, "%f", 1);
  dQhr(ii,cycle) = fscanf(fid, "%f", 1);

  # Reads line 3
  xb(ii,cycle) = fscanf(fid, "%f", 1);
  xbdot(ii,cycle) = fscanf(fid, "%f", 1);

  # Reads line 4
  m_in(ii,cycle)    = fscanf(fid, "%f", 1);
  m_out(ii,cycle)   = fscanf(fid, "%f", 1);
  vol(ii,cycle)     = fscanf(fid, "%f", 1);
  vdot(ii,cycle)    = fscanf(fid, "%f", 1);
  m_fuel(ii,cycle)  = fscanf(fid, "%f", 1);
  m_air(ii,cycle)   = fscanf(fid, "%f", 1);
  m_res(ii,cycle)   = fscanf(fid, "%f", 1);
  dQ_ht(ii,cycle)   = fscanf(fid, "%f", 1);
  dQ_chem(ii,cycle) = fscanf(fid, "%f", 1);
  bas = fscanf(fid, "%s", 1);

end
fclose(fid);
