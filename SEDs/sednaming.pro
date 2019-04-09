pro sednaming

ff='Tseds.dat'
readcol,ff,teff,logg,zzo

nn=(size(teff))[1]

for ii=0,nn-1 do begin
  if (logg[ii] eq fix(logg[ii]) ) then begin
    print,round(teff[ii]),logg[ii],round(zzo[ii]),'   TlustyOB/t'+strcompress(round(teff[ii]),/remove_all)+'.g'+strcompress(fix(logg[ii]),/remove_all)+'.'+strcompress(fix(100*(logg[ii]-fix(logg[ii]))),/remove_all)+'0z1.dat.txt'
  endif else print,round(teff[ii]),logg[ii],round(zzo[ii]),'   TlustyOB/t'+strcompress(round(teff[ii]),/remove_all)+'.g'+strcompress(fix(logg[ii]),/remove_all)+'.'+strcompress(fix(100*(logg[ii]-fix(logg[ii]))),/remove_all)+'z1.dat.txt'
endfor



stop
end