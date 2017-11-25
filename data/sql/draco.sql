SELECT
 p.ra, p.dec, p.psfMag_g, psfMagErr_g, p.psfMag_r, p.psfMagErr_r, p.psfMag_i, p.psfMagErr_i
FROM PhotoObjAll AS p
WHERE
  p.ra BETWEEN 259.05 AND 261.05
  AND p.dec BETWEEN 56.92 AND 58.92
  AND dbo.fDistanceArcMinEq(260.05, 57.92, p.ra, p.dec) < 30.
  AND p.clean = 1
