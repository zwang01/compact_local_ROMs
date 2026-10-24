function rom = eff_enrich(rom, S, par, U, W, dvec)
% Reuse the existing POD/window policy and cache exact convection stencils
% applied to the new state basis. All preparation is timed by the driver.
rom=brinkman_rom_enrich(rom,S,U,W);
if strcmp(rom.eff_variant,'baseline'), return; end
if ~isfield(rom,'C'), rom.C=eff_convection_maps(par,S); end
C=rom.C;
rom.X=C.X*rom.Phi; rom.Y=C.Y*rom.Phi;
rom.Dx=C.Dx*rom.Phi; rom.Dy=C.Dy*rom.Phi;
if strcmp(rom.eff_variant,'sampled'), rom=eff_sample_prepare(rom,S,dvec,U); end
end
