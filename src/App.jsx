import React, { useState, useEffect, useCallback } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, ReferenceLine, Area, LabelList, Customized } from 'recharts';

function App() {
  // State untuk domain YAxis SFD dan BMD
  const [sfdYAxisDomain, setSfdYAxisDomain] = useState(['auto', 'auto']);
  const [bmdYAxisDomain, setBmdYAxisDomain] = useState(['auto', 'auto']);

  // State untuk ticks sumbu Y diagram
  const [sfdYAxisTicks, setSfdYAxisTicks] = useState([]);
  const [bmdYAxisTicks, setBmdYAxisTicks] = useState([]);

  // State untuk data balok dan beban
  const [beamLength, setBeamLength] = useState(10); // Panjang balok dalam meter
  const [concentratedLoads, setConcentratedLoads] = useState([]); // Array beban terpusat { magnitude, position, unit }
  const [distributedLoads, setDistributedLoads] = useState([]); // Array beban terbagi merata { intensity, start, end, unit }

  // State untuk tumpuan: selalu ada dua tumpuan
  const [supports, setSupports] = useState([
    { id: 'support-A', type: 'hinge', position: 0 },
    { id: 'support-B', type: 'roller', position: 10 } // Initial position, will be updated by useEffect
  ]);

  // State untuk hasil perhitungan
  const [reactions, setReactions] = useState({}); // Gaya reaksi tumpuan, disimpan berdasarkan ID tumpuan
  const [sfdData, setSfdData] = useState([]); // Data untuk diagram gaya geser (lebih detail untuk plotting)
  const [bmdData, setBmdData] = useState([]); // Data untuk diagram momen lentur (lebih detail untuk plotting)
  const [tableSfdData, setTableSfdData] = useState([]); // Data untuk tabel gaya geser (interval 1m)
  const [tableBmdData, setTableBmdData] = useState([]); // Data untuk tabel momen lentur (interval 1m)
  const [errorMessage, setErrorMessage] = useState(''); // Pesan error

  // State untuk unit output yang dipilih
  const [outputUnit, setOutputUnit] = useState('kN');

  // State untuk label sumbu Y diagram
  const [sfdYAxisLabel, setSfdYAxisLabel] = useState('Gaya Geser (kN)');
  const [bmdYAxisLabel, setBmdYAxisLabel] = useState('Momen (kNm)');


  // Konversi unit (faktor untuk mengkonversi unit ke kN)
  const unitConversionFactor = {
    'kN': 1,
    'ton': 9.80665, // 1 metric ton = 9.80665 kN
    'kg': 0.00980665 // 1 kg = 0.00980665 kN
  };

  // Efek samping untuk mengupdate posisi tumpuan B jika panjang balok berubah
  useEffect(() => {
    setSupports(prevSupports => {
      const updatedSupports = prevSupports.map(s =>
        s.id === 'support-B' ? { ...s, position: beamLength } : s
      );
      // Sort supports by position for consistent display and calculation order
      return updatedSupports.sort((a, b) => a.position - b.position);
    });
  }, [beamLength]);

  // --- Fungsi Penambahan/Penghapusan Beban ---

  const addConcentratedLoad = () => {
    setConcentratedLoads([...concentratedLoads, { magnitude: 10, position: beamLength / 2, unit: 'kN' }]);
  };

  const removeConcentratedLoad = (index) => {
    const newLoads = concentratedLoads.filter((_, i) => i !== index);
    setConcentratedLoads(newLoads);
  };

  const updateConcentratedLoad = (index, field, value) => {
    const newLoads = [...concentratedLoads];
    newLoads[index][field] = field === 'magnitude' || field === 'position' ? parseFloat(value) : value;
    setConcentratedLoads(newLoads);
  };

  const addDistributedLoad = () => {
    setDistributedLoads([...distributedLoads, { intensity: 5, start: beamLength / 4, end: beamLength * 3 / 4, unit: 'kN/m' }]);
  };

  const removeDistributedLoad = (index) => {
    const newLoads = distributedLoads.filter((_, i) => i !== index);
    setDistributedLoads(newLoads);
  };

  const updateDistributedLoad = (index, field, value) => {
    const newLoads = [...distributedLoads];
    newLoads[index][field] = field === 'intensity' || field === 'start' || field === 'end' ? parseFloat(value) : value;
    setDistributedLoads(newLoads);
  };

  // --- Fungsi Update Tumpuan ---
  const updateSupport = (id, field, value) => {
    setSupports(prevSupports => {
      let newSupports = prevSupports.map(s =>
        s.id === id ? { ...s, [field]: field === 'position' ? parseFloat(value) : value } : s
      );

      // Enforce one hinge and one roller for static determinacy
      if (field === 'type') {
        const targetSupport = newSupports.find(s => s.id === id);
        const otherSupport = newSupports.find(s => s.id !== id);

        // If changing type makes both supports the same type, flip the other one
        if (targetSupport.type === otherSupport.type) {
          otherSupport.type = (targetSupport.type === 'hinge' ? 'roller' : 'hinge');
        }
      }

      // Sort supports by position for consistent A/B identification
      return newSupports.sort((a, b) => a.position - b.position);
    });
  };

  // --- Fungsi untuk menghasilkan nilai ticks sumbu Y yang rapi ---
  const generateNiceTicksForYAxis = (minVal, maxVal) => {
    const ticks = [];
    // Handle flat line case (minVal and maxVal are very close or equal)
    if (Math.abs(maxVal - minVal) < 1e-6) { // Use epsilon for floating point comparison
      const singleTick = parseFloat(minVal.toFixed(2));
      ticks.push(singleTick);
      if (singleTick !== 0) ticks.push(0); // Ensure 0 is included if not flat at 0
      return [...new Set(ticks)].sort((a, b) => a - b);
    }

    const paddingRatio = 0.1; // 10% padding on each side
    const paddedMin = minVal - (maxVal - minVal) * paddingRatio;
    const paddedMax = maxVal + (maxVal - minVal) * paddingRatio;
    const span = paddedMax - paddedMin;

    const numDesiredTicks = 5; // Target jumlah ticks yang diinginkan

    // Tentukan ukuran langkah yang 'bagus'
    let roughStep = span / numDesiredTicks;
    let exponent = Math.floor(Math.log10(roughStep));
    let fraction = roughStep / Math.pow(10, exponent);
    let niceFraction;

    if (fraction <= 1) niceFraction = 1;
    else if (fraction <= 2) niceFraction = 2;
    else if (fraction <= 5) niceFraction = 5;
    else niceFraction = 10;

    const step = niceFraction * Math.pow(10, exponent);

    // Mulai dari tick pertama yang kurang dari atau sama dengan minVal, tetapi disesuaikan
    let currentTick = Math.floor(paddedMin / step) * step;

    const epsilon = 1e-9; // Use a smaller epsilon for more robust float comparison

    // Tambahkan tick hingga melewati maxVal
    while (currentTick <= paddedMax + epsilon) { // Loop until we pass the padded max
      // Pastikan tick diformat dengan presisi yang wajar untuk menghindari masalah floating point
      const formattedTick = parseFloat(currentTick.toFixed(Math.max(0, -exponent) + 2)); // Adjust precision based on step
      ticks.push(formattedTick);
      currentTick += step;
    }

    // Pastikan 0 disertakan secara tepat jika berada dalam rentang
    if (minVal < 0 && maxVal > 0) {
      let foundZero = false;
      for(const tick of ticks) {
          if (Math.abs(tick) < 1e-9) { // Check if a tick is very close to zero
              foundZero = true;
              break;
          }
      }
      if (!foundZero) ticks.push(0);
    }

    // Urutkan dan hapus duplikat (karena epsilon atau penyisipan 0)
    return [...new Set(ticks)].sort((a, b) => a - b);
  };


  // --- Fungsi Perhitungan Utama ---
  const calculateBeam = useCallback(() => {
    // Validasi input balok
    if (beamLength <= 0 || isNaN(beamLength)) {
      setErrorMessage('Panjang balok harus lebih besar dari 0.');
      setReactions({});
      setSfdData([]);
      setBmdData([]);
      setTableSfdData([]);
      setTableBmdData([]);
      return;
    }

    // Validasi posisi tumpuan dan jenis
    const [supportA, supportB] = supports.sort((a, b) => a.position - b.position); // Ensure A is always left of B

    if (isNaN(supportA.position) || isNaN(supportB.position) || supportA.position < 0 || supportB.position > beamLength || supportA.position >= supportB.position) {
      setErrorMessage('Posisi tumpuan tidak valid. Pastikan mereka dalam rentang balok dan tumpuan kiri berada di sebelah kiri tumpuan kanan.');
      setReactions({});
      setSfdData([]);
      setBmdData([]);
      setTableSfdData([]);
      setTableBmdData([]);
      return;
    }

    if (supportA.type === supportB.type) {
        setErrorMessage('Konfigurasi tumpuan tidak valid. Harus ada tepat satu tumpuan sendi dan satu tumpuan rol.');
        setReactions({});
        setSfdData([]);
        setBmdData([]);
        setTableSfdData([]);
        setTableBmdData([]);
        return;
    }

    // Validasi beban terpusat
    for (const load of concentratedLoads) {
      if (isNaN(load.magnitude) || isNaN(load.position) || load.position < 0 || load.position > beamLength || !unitConversionFactor[load.unit]) {
        setErrorMessage('Beban terpusat tidak valid (nilai/posisi/unit). Pastikan posisi dalam rentang balok.');
        setReactions({});
        setSfdData([]);
        setBmdData([]);
        setTableSfdData([]);
        setTableBmdData([]);
        return;
      }
    }

    // Validasi beban terbagi merata
    for (const load of distributedLoads) {
      if (isNaN(load.intensity) || isNaN(load.start) || isNaN(load.end) || load.start < 0 || load.end > beamLength || load.start >= load.end || !unitConversionFactor[load.unit.split('/')[0]]) {
        setErrorMessage('Beban terbagi merata tidak valid (nilai/posisi/unit). Pastikan posisi dalam rentang balok dan titik awal lebih kecil dari titik akhir.');
        setReactions({});
        setSfdData([]);
        setBmdData([]);
        setTableSfdData([]);
        setTableBmdData([]);
        return;
      }
    }

    setErrorMessage(''); // Hapus pesan error jika validasi berhasil

    // 1. Hitung Gaya Reaksi Tumpuan (semua dalam kN)
    let hingeSupportCalc = supports.find(s => s.type === 'hinge');
    let rollerSupportCalc = supports.find(s => s.type === 'roller');

    // Default ke supportA dan supportB jika jenis tidak ditemukan atau tidak ada 2 tumpuan
    if (!hingeSupportCalc && supports.length > 0) hingeSupportCalc = supports[0];
    if (!rollerSupportCalc && supports.length > 1) rollerSupportCalc = supports[1];

    // Fallback jika masih null (misal, hanya ada 1 support atau tidak ada)
    if (!hingeSupportCalc) hingeSupportCalc = supports[0];
    if (!rollerSupportCalc) rollerSupportCalc = supports[1];


    const spanLength = rollerSupportCalc.position - hingeSupportCalc.position;

    let sumMomentsAboutHinge = 0;

    // Momen akibat beban terpusat terhadap tumpuan sendi
    concentratedLoads.forEach(load => {
      const magnitudeInKN = load.magnitude * unitConversionFactor[load.unit];
      sumMomentsAboutHinge += magnitudeInKN * (load.position - hingeSupportCalc.position);
    });

    // Momen akibat beban terbagi merata terhadap tumpuan sendi
    distributedLoads.forEach(load => {
      const intensityInKN_m = load.intensity * unitConversionFactor[load.unit.split('/')[0]];
      // Perhatikan hanya bagian beban yang berada dalam rentang balok
      const loadStartInSpan = Math.max(load.start, 0);
      const loadEndInSpan = Math.min(load.end, beamLength);

      if (loadEndInSpan > loadStartInSpan) {
        // Hitung momen dari beban terbagi merata
        // Potongan UDL yang aktif dihitung relatif terhadap tumpuan sendi
        const effectiveStartRelHinge = Math.max(0, loadStartInSpan - hingeSupportCalc.position);
        const effectiveEndRelHinge = Math.min(spanLength, loadEndInSpan - hingeSupportCalc.position);

        if (effectiveEndRelHinge > effectiveStartRelHinge) {
          const effectiveLoadLength = effectiveEndRelHinge - effectiveStartRelHinge;
          const equivalentForce = intensityInKN_m * effectiveLoadLength;
          const centroidPositionRelHinge = effectiveStartRelHinge + effectiveLoadLength / 2;
          sumMomentsAboutHinge += equivalentForce * centroidPositionRelHinge;
        }
      }
    });

    const R_roller_V_inKN = sumMomentsAboutHinge / spanLength;

    let totalExternalForceInKN = 0;
    concentratedLoads.forEach(load => { totalExternalForceInKN += load.magnitude * unitConversionFactor[load.unit]; });
    distributedLoads.forEach(load => {
        const loadLength = load.end - load.start;
        totalExternalForceInKN += load.intensity * unitConversionFactor[load.unit.split('/')[0]] * loadLength;
    });

    const R_hinge_V_inKN = totalExternalForceInKN - R_roller_V_inKN;

    // Konversi hasil reaksi ke outputUnit
    const displayReactions = {
      [hingeSupportCalc.id]: (R_hinge_V_inKN / unitConversionFactor[outputUnit]).toFixed(2),
      [rollerSupportCalc.id]: (R_roller_V_inKN / unitConversionFactor[outputUnit]).toFixed(2)
    };
    setReactions(displayReactions);

    // Helper functions for shear and moment calculation in KN and KN-m
    const calculateShearForceAtX = (x_coord) => {
      let shear = 0;
      const calcEpsilon = 0.0001; // Epsilon for calculation points

      // Add hinge reaction
      if (x_coord >= hingeSupportCalc.position - calcEpsilon) {
        shear += R_hinge_V_inKN;
      }
      // Add roller reaction
      if (x_coord >= rollerSupportCalc.position - calcEpsilon) {
        shear += R_roller_V_inKN;
      }

      // Subtract concentrated loads
      concentratedLoads.forEach(load => {
        if (x_coord >= load.position - calcEpsilon) {
          shear -= load.magnitude * unitConversionFactor[load.unit];
        }
      });

      // Subtract distributed loads
      distributedLoads.forEach(load => {
        const intensityInKN_m = load.intensity * unitConversionFactor[load.unit.split('/')[0]];
        const activeStart = Math.max(0, load.start);
        const activeEnd = Math.min(beamLength, load.end);

        if (x_coord > activeStart - calcEpsilon) {
          const activeLength = Math.min(x_coord, activeEnd) - activeStart;
          if (activeLength > 0) {
            shear -= intensityInKN_m * activeLength;
          }
        }
      });
      return shear;
    };

    const calculateMomentAtX = (x_coord) => {
        let moment = 0;
        const calcEpsilon = 0.0001;

        // Moment due to hinge reaction
        if (x_coord >= hingeSupportCalc.position - calcEpsilon) {
            moment += R_hinge_V_inKN * (x_coord - hingeSupportCalc.position);
        }
        // Moment due to roller reaction
        if (x_coord >= rollerSupportCalc.position - calcEpsilon) {
            moment += R_roller_V_inKN * (x_coord - rollerSupportCalc.position);
        }

        // Moment due to concentrated loads
        concentratedLoads.forEach(load => {
            if (x_coord >= load.position - calcEpsilon) {
                moment -= load.magnitude * unitConversionFactor[load.unit] * (x_coord - load.position);
            }
        });

        // Moment due to distributed loads
        distributedLoads.forEach(load => {
            const intensityInKN_m = load.intensity * unitConversionFactor[load.unit.split('/')[0]];
            const activeStart = Math.max(0, load.start);
            const activeEnd = Math.min(beamLength, load.end);

            if (x_coord > activeStart - calcEpsilon) {
                const activeLength = Math.min(x_coord, activeEnd) - activeStart;
                if (activeLength > 0) {
                    moment -= (intensityInKN_m * activeLength) * (x_coord - (activeStart + activeLength / 2));
                }
            }
        });
        return moment;
    };

    const epsilonForPlotting = 0.0001; // Epsilon untuk menangani diskontinuitas visual pada grafik

    // 1. Generate Data for SFD GRAPH (detailed plotting points including epsilon for jumps)
    // Kumpulkan semua titik kritis (posisi beban, posisi tumpuan, 0, beamLength)
    const allSignificantXCoords = new Set([0, beamLength]);
    supports.forEach(s => allSignificantXCoords.add(s.position));
    concentratedLoads.forEach(l => allSignificantXCoords.add(l.position));
    distributedLoads.forEach(l => {
      allSignificantXCoords.add(l.start);
      allSignificantXCoords.add(l.end);
    });

    // Tambahkan titik sedikit sebelum dan sesudah titik kritis untuk menangkap diskontinuitas
    const criticalPlotPoints = Array.from(allSignificantXCoords).sort((a, b) => a - b);
    const plotXPoints = new Set();

    criticalPlotPoints.forEach(point => {
      plotXPoints.add(point);
      // Only add epsilon points if they are within the beam length and distinct
      if (point > 0) plotXPoints.add(Math.max(0, point - epsilonForPlotting));
      if (point < beamLength) plotXPoints.add(Math.min(beamLength, point + epsilonForPlotting));
    });

    // Tambahkan titik menengah untuk kurva yang halus
    const sortedUniquePlotPoints = Array.from(plotXPoints).sort((a, b) => a - b);
    const finalPlottingPoints = new Set();

    const numIntermediateSegments = 50; // Jumlah segmen untuk setiap interval kritis
    sortedUniquePlotPoints.forEach((point, index) => {
      finalPlottingPoints.add(point);
      if (index < sortedUniquePlotPoints.length - 1) {
        const nextPoint = sortedUniquePlotPoints[index + 1];
        const segmentLength = nextPoint - point;
        if (segmentLength > 0) { // Avoid division by zero for identical points
            const stepSize = segmentLength / numIntermediateSegments;
            for (let i = 1; i < numIntermediateSegments; i++) {
              finalPlottingPoints.add(point + (i * stepSize));
            }
        }
      }
    });
    finalPlottingPoints.add(beamLength); // Ensure beamLength is included

    const sortedFinalPlottingPoints = Array.from(finalPlottingPoints).sort((a, b) => a - b);

    const tempSfdDataPlot = []; // Data untuk plotting SFD
    const tempBmdDataPlot = []; // Data untuk plotting BMD

    let minShearPlot_inKN = 0, maxShearPlot_inKN = 0;
    let minMomentPlot_inKNm = 0, maxMomentPlot_inKNm = 0;

    sortedFinalPlottingPoints.forEach(x => {
      let currentShear_inKN = calculateShearForceAtX(x);
      let currentMoment_inKNm = calculateMomentAtX(x);

      // Special handling for the very last point on SFD graph:
      // If the last point is at beamLength and has a support,
      // we want the shear value *just before* that support reaction is applied.
      const supportAtBeamEnd = supports.find(s => Math.abs(s.position - beamLength) < epsilonForPlotting);
      if (Math.abs(x - beamLength) < epsilonForPlotting && supportAtBeamEnd) {
           currentShear_inKN = calculateShearForceAtX(beamLength - epsilonForPlotting);
           // Ensure it's not a tiny negative due to floating point if it should be zero
           if (Math.abs(currentShear_inKN) < 1e-9) currentShear_inKN = 0;
      }

      const shearValueDisplay = parseFloat((currentShear_inKN / unitConversionFactor[outputUnit]).toFixed(2));
      const momentValueDisplay = parseFloat((currentMoment_inKNm / unitConversionFactor[outputUnit]).toFixed(2));

      tempSfdDataPlot.push({
        x: parseFloat(x.toFixed(2)),
        V: shearValueDisplay,
        V_pos: shearValueDisplay >= 0 ? shearValueDisplay : null, // For positive area fill
        V_neg: shearValueDisplay <= 0 ? shearValueDisplay : null  // For negative area fill
      });
      tempBmdDataPlot.push({ x: parseFloat(x.toFixed(2)), M: momentValueDisplay });

      minShearPlot_inKN = Math.min(minShearPlot_inKN, currentShear_inKN);
      maxShearPlot_inKN = Math.max(maxShearPlot_inKN, currentShear_inKN);
      minMomentPlot_inKNm = Math.min(minMomentPlot_inKNm, currentMoment_inKNm);
      maxMomentPlot_inKNm = Math.max(maxMomentPlot_inKNm, currentMoment_inKNm);
    });

    // Set data untuk grafik
    setSfdData(tempSfdDataPlot); // Menggunakan data detail untuk grafik SFD
    setBmdData(tempBmdDataPlot); // Menggunakan data detail untuk grafik BMD


    // 2. Generate Data for TABLES (integer intervals with specific end point handling for SFD)
    const newTableSfdData = [];
    const newTableBmdData = [];
    const epsilonForTableCheck = 0.0001; // Epsilon khusus untuk pengecekan tabel

    for (let i = 0; i <= beamLength; i++) {
        let x = i;
        let currentShear_inKN_table = calculateShearForceAtX(x); // Hitung awal termasuk semua gaya di X

        // Penanganan khusus untuk titik terakhir di tabel SFD:
        // Kita ingin nilai gaya geser *sebelum* gaya reaksi di beamLength bekerja
        const supportAtBeamEndTable = supports.find(s => Math.abs(s.position - beamLength) < epsilonForTableCheck);

        if (Math.abs(x - beamLength) < epsilonForTableCheck && supportAtBeamEndTable) {
            let reactionToSubtract = 0;
            // Tentukan reaksi tumpuan mana yang akan dikurangi (sendi atau rol)
            if (supportAtBeamEndTable.id === hingeSupportCalc.id) {
                reactionToSubtract = R_hinge_V_inKN;
            } else if (supportAtBeamEndTable.id === rollerSupportCalc.id) {
                reactionToSubtract = R_roller_V_inKN;
            }
            currentShear_inKN_table -= reactionToSubtract;

            // Pastikan pembulatan ke nol jika nilai sangat dekat nol
            if (Math.abs(currentShear_inKN_table) < 1e-9) currentShear_inKN_table = 0;
        }

        let currentMoment_inKNm_table = calculateMomentAtX(x); // Momen tetap dihitung pada X

        newTableSfdData.push({
            x: i, // Posisi integer untuk tabel
            V: parseFloat((currentShear_inKN_table / unitConversionFactor[outputUnit]).toFixed(2))
        });
        newTableBmdData.push({
            x: i, // Posisi integer untuk tabel
            M: parseFloat((currentMoment_inKNm_table / unitConversionFactor[outputUnit]).toFixed(2))
        });
    }
    setTableSfdData(newTableSfdData);
    setTableBmdData(newTableBmdData);


    // 3. Set Y-Axis Domains and Ticks for SFD (based on detailed graph data)
    const minShearDisplayForGraph = minShearPlot_inKN / unitConversionFactor[outputUnit];
    const maxShearDisplayForGraph = maxShearPlot_inKN / unitConversionFactor[outputUnit];
    let sfdTicks = generateNiceTicksForYAxis(minShearDisplayForGraph, maxShearDisplayForGraph);
    setSfdYAxisTicks(sfdTicks);
    setSfdYAxisDomain(sfdTicks.length > 0 ? [sfdTicks[0], sfdTicks[sfdTicks.length - 1]] : ['auto', 'auto']);

    // 4. Set Y-Axis Domains and Ticks for BMD (based on detailed plot data)
    const minMomentDisplayDetailed = minMomentPlot_inKNm / unitConversionFactor[outputUnit];
    const maxMomentDisplayDetailed = maxMomentPlot_inKNm / unitConversionFactor[outputUnit];
    let bmdTicks = generateNiceTicksForYAxis(minMomentDisplayDetailed, maxMomentDisplayDetailed);
    setBmdYAxisTicks(bmdTicks);
    setBmdYAxisDomain(bmdTicks.length > 0 ? [bmdTicks[0], bmdTicks[bmdTicks.length - 1]] : ['auto', 'auto']);

    // Update label sumbu Y
    setSfdYAxisLabel(`Gaya Geser (${outputUnit})`);
    const momentUnitLabel = outputUnit === 'kN' ? 'kNm' : `${outputUnit}m`;
    setBmdYAxisLabel(`Momen (${momentUnitLabel})`);


  }, [beamLength, concentratedLoads, distributedLoads, supports, outputUnit]); // Removed 'reactions' from dependencies

  // X-axis ticks generation (dari permintaan sebelumnya)
  const generateXAxisTicks = useCallback(() => {
    const ticks = [];
    for (let i = 0; i <= beamLength; i++) {
      ticks.push(i);
    }
    return ticks;
  }, [beamLength]);

  // Efek samping untuk menghitung ulang saat data atau tumpuan berubah
  useEffect(() => {
    calculateBeam();
  }, [calculateBeam]);

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100 p-4 font-inter text-gray-800">
      <div className="container mx-auto bg-white shadow-xl rounded-2xl p-8 border border-blue-200">
        <h1 className="text-4xl font-extrabold text-center text-blue-700 mb-2 tracking-tight">
          Kalkulator Analisis Balok Sederhana
        </h1>
        <p className="text-center text-red-600 text-sm mb-8">Dibuat oleh ZUELCAR'9</p>

        {errorMessage && (
          <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded-lg relative mb-6" role="alert">
            <strong className="font-bold">Error!</strong>
            <span className="block sm:inline ml-2">{errorMessage}</span>
          </div>
        )}

        {/* Bagian Input */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-8 mb-10">
          {/* Input Panjang Balok */}
          <div className="p-6 bg-blue-50 rounded-xl shadow-md border border-blue-200">
            <h2 className="text-2xl font-semibold text-blue-600 mb-4">Properti Balok</h2>
            <div className="mb-4">
              <label htmlFor="beamLength" className="block text-sm font-medium text-gray-700 mb-1">
                Panjang Balok (L):
              </label>
              <div className="relative mt-1 rounded-md shadow-sm">
                <input
                  type="number"
                  id="beamLength"
                  value={beamLength}
                  onChange={(e) => setBeamLength(parseFloat(e.target.value))}
                  min="0.1"
                  step="0.1"
                  className="block w-full rounded-lg border-gray-300 pr-10 focus:border-blue-500 focus:ring-blue-500 shadow-sm sm:text-sm p-2.5"
                  placeholder="Masukkan panjang"
                />
                <div className="pointer-events-none absolute inset-y-0 right-0 flex items-center pr-3">
                  <span className="text-gray-500 sm:text-sm">m</span>
                </div>
              </div>
            </div>
          </div>

          {/* Input Tumpuan */}
          <div className="p-6 bg-yellow-50 rounded-xl shadow-md border border-yellow-200">
            <h2 className="text-2xl font-semibold text-yellow-700 mb-4">Tumpuan</h2>
            {supports.map((s, index) => (
              <div key={s.id} className="flex items-center space-x-3 mb-3 bg-yellow-100 p-3 rounded-lg shadow-sm">
                <div className="flex-grow grid grid-cols-2 gap-3">
                  <div>
                    <label htmlFor={`support-type-${s.id}`} className="block text-xs font-medium text-gray-700 mb-0.5">
                      Jenis Tumpuan ({s.id.split('-')[1]}):
                    </label>
                    <select
                      id={`support-type-${s.id}`}
                      value={s.type}
                      onChange={(e) => updateSupport(s.id, 'type', e.target.value)}
                      className="block w-full rounded-md border-gray-300 focus:border-yellow-500 focus:ring-yellow-500 sm:text-sm p-2"
                    >
                      <option value="hinge">Sendi</option>
                      <option value="roller">Rol</option>
                    </select>
                  </div>
                  <div>
                    <label htmlFor={`support-pos-${s.id}`} className="block text-xs font-medium text-gray-700 mb-0.5">
                      Posisi dari Kiri (m):
                    </label>
                    <input
                      type="number"
                      id={`support-pos-${s.id}`}
                      value={s.position}
                      onChange={(e) => updateSupport(s.id, 'position', e.target.value)}
                      min="0"
                      max={beamLength}
                      step="0.1"
                      className="block w-full rounded-md border-gray-300 focus:border-yellow-500 focus:ring-yellow-500 sm:text-sm p-2"
                    />
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>

        {/* Input Beban Terpusat */}
        <div className="p-6 bg-green-50 rounded-xl shadow-md border border-green-200 mb-10">
          <h2 className="text-2xl font-semibold text-green-600 mb-4">Beban Terpusat</h2>
          {concentratedLoads.map((load, index) => (
            <div key={index} className="flex items-center space-x-3 mb-3 bg-green-100 p-3 rounded-lg shadow-sm">
              <div className="flex-grow grid grid-cols-3 gap-3">
                <div>
                  <label htmlFor={`P-magnitude-${index}`} className="block text-xs font-medium text-gray-700 mb-0.5">
                    Besar:
                  </label>
                  <input
                    type="number"
                    id={`P-magnitude-${index}`}
                    value={load.magnitude}
                    onChange={(e) => updateConcentratedLoad(index, 'magnitude', e.target.value)}
                    step="0.1"
                    className="block w-full rounded-md border-gray-300 focus:border-green-500 focus:ring-green-500 sm:text-sm p-2"
                  />
                </div>
                <div>
                  <label htmlFor={`P-unit-${index}`} className="block text-xs font-medium text-gray-700 mb-0.5">
                    Unit:
                  </label>
                  <select
                    id={`P-unit-${index}`}
                    value={load.unit}
                    onChange={(e) => updateConcentratedLoad(index, 'unit', e.target.value)}
                    className="block w-full rounded-md border-gray-300 focus:border-green-500 focus:ring-green-500 sm:text-sm p-2"
                  >
                    <option value="kN">kN</option>
                    <option value="ton">ton</option>
                    <option value="kg">kg</option>
                  </select>
                </div>
                <div>
                  <label htmlFor={`P-position-${index}`} className="block text-xs font-medium text-gray-700 mb-0.5">
                    Posisi dari Kiri (m):
                  </label>
                  <input
                    type="number"
                    id={`P-position-${index}`}
                    value={load.position}
                    onChange={(e) => updateConcentratedLoad(index, 'position', e.target.value)}
                    min="0"
                    max={beamLength}
                    step="0.1"
                    className="block w-full rounded-md border-gray-300 focus:border-green-500 focus:ring-green-500 sm:text-sm p-2"
                  />
                </div>
              </div>
              <button
                onClick={() => removeConcentratedLoad(index)}
                className="bg-red-500 hover:bg-red-600 text-white p-2 rounded-full shadow-md transition-all duration-200 ease-in-out"
              >
                <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
                  <path fillRule="evenodd" d="M9 2a1 1 0 00-.894.553L7.382 4H4a1 1 0 000 2v10a2 2 0 002 2h8a2 2 0 002-2V6a1 1 0 100-2h-3.382l-.724-1.447A1 1 0 0011 2H9zM7 8a1 1 0 012 0v6a1 1 0 11-2 0V8zm6 0a1 1 0 11-2 0v6a1 1 0 112 0V8z" clipRule="evenodd" />
                </svg>
              </button>
            </div>
          ))}
          <button
            onClick={addConcentratedLoad}
            className="mt-4 w-full bg-green-600 hover:bg-green-700 text-white font-bold py-2 px-4 rounded-lg shadow-md transition-all duration-200 ease-in-out"
          >
            Tambah Beban Terpusat
          </button>
        </div>

        {/* Input Beban Terbagi Merata */}
        <div className="p-6 bg-purple-50 rounded-xl shadow-md border border-purple-200 mb-10">
          <h2 className="text-2xl font-semibold text-purple-600 mb-4">Beban Terbagi Merata</h2>
          {distributedLoads.map((load, index) => (
            <div key={index} className="flex items-center space-x-3 mb-3 bg-purple-100 p-3 rounded-lg shadow-sm">
              <div className="flex-grow grid grid-cols-4 gap-3">
                <div>
                  <label htmlFor={`UDL-intensity-${index}`} className="block text-xs font-medium text-gray-700 mb-0.5">
                    Intensitas:
                  </label>
                  <input
                    type="number"
                    id={`UDL-intensity-${index}`}
                    value={load.intensity}
                    onChange={(e) => updateDistributedLoad(index, 'intensity', e.target.value)}
                    step="0.1"
                    className="block w-full rounded-md border-gray-300 focus:border-purple-500 focus:ring-purple-500 sm:text-sm p-2"
                  />
                </div>
                <div>
                  <label htmlFor={`UDL-unit-${index}`} className="block text-xs font-medium text-gray-700 mb-0.5">
                    Unit:
                  </label>
                  <select
                    id={`UDL-unit-${index}`}
                    value={load.unit}
                    onChange={(e) => updateDistributedLoad(index, 'unit', e.target.value)}
                    className="block w-full rounded-md border-gray-300 focus:border-purple-500 focus:ring-purple-500 sm:text-sm p-2"
                  >
                    <option value="kN/m">kN/m</option>
                    <option value="ton/m">ton/m</option>
                    <option value="kg/m">kg/m</option>
                  </select>
                </div>
                <div>
                  <label htmlFor={`UDL-start-${index}`} className="block text-xs font-medium text-gray-700 mb-0.5">
                    Mulai dari (m):
                  </label>
                  <input
                    type="number"
                    id={`UDL-start-${index}`}
                    value={load.start}
                    onChange={(e) => updateDistributedLoad(index, 'start', e.target.value)}
                    min="0"
                    max={beamLength}
                    step="0.1"
                    className="block w-full rounded-md border-gray-300 focus:border-purple-500 focus:ring-purple-500 sm:text-sm p-2"
                  />
                </div>
                <div>
                  <label htmlFor={`UDL-end-${index}`} className="block text-xs font-medium text-gray-700 mb-0.5">
                    Berakhir di (m):
                  </label>
                  <input
                    type="number"
                    id={`UDL-end-${index}`}
                    value={load.end}
                    onChange={(e) => updateDistributedLoad(index, 'end', e.target.value)}
                    min="0"
                    max={beamLength}
                    step="0.1"
                    className="block w-full rounded-md border-gray-300 focus:border-purple-500 focus:ring-purple-500 sm:text-sm p-2"
                  />
                </div>
              </div>
              <button
                onClick={() => removeDistributedLoad(index)}
                className="bg-red-500 hover:bg-red-600 text-white p-2 rounded-full shadow-md transition-all duration-200 ease-in-out"
              >
                <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
                  <path fillRule="evenodd" d="M9 2a1 1 0 00-.894.553L7.382 4H4a1 1 0 000 2v10a2 2 0 002 2h8a2 2 0 002-2V6a1 1 0 100-2h-3.382l-.724-1.447A1 1 0 0011 2H9zM7 8a1 1 0 012 0v6a1 1 0 11-2 0V8zm6 0a1 1 0 11-2 0v6a1 1 0 112 0V8z" clipRule="evenodd" />
                </svg>
              </button>
            </div>
          ))}
          <button
            onClick={addDistributedLoad}
            className="mt-4 w-full bg-purple-600 hover:bg-purple-700 text-white font-bold py-2 px-4 rounded-lg shadow-md transition-all duration-200 ease-in-out"
          >
            Tambah Beban Terbagi Merata
          </button>
        </div>

        {/* Tombol Hitung */}
        <div className="mb-10 text-center">
          <button
            onClick={calculateBeam}
            className="bg-blue-600 hover:bg-blue-700 text-white font-bold py-3 px-8 rounded-full shadow-lg transition-all duration-300 ease-in-out transform hover:scale-105"
          >
            Hitung Gaya Reaksi & Diagram
          </button>
        </div>

        {/* Hasil Perhitungan */}
        <div className="p-6 bg-yellow-50 rounded-xl shadow-md border border-yellow-200 mb-10">
          <h2 className="text-2xl font-semibold text-yellow-700 mb-4 text-center">Hasil Perhitungan</h2>
          <div className="mb-4">
            <label htmlFor="outputUnit" className="block text-sm font-medium text-gray-700 mb-1">
              Tampilkan Hasil dalam Unit:
            </label>
            <select
              id="outputUnit"
              value={outputUnit}
              onChange={(e) => setOutputUnit(e.target.value)}
              className="block w-full rounded-md border-gray-300 focus:border-blue-500 focus:ring-blue-500 sm:text-sm p-2"
            >
              <option value="kN">kN</option>
              <option value="ton">ton</option>
              <option value="kg">kg</option>
            </select>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-lg font-medium">
            {supports.map(s => (
              <div key={s.id} className="bg-yellow-100 p-4 rounded-lg shadow-sm text-center">
                <p>Gaya Reaksi {s.id.split('-')[1]} ({s.type === 'hinge' ? 'Sendi' : 'Rol'})
                  <span className="font-bold text-yellow-800">
                    {reactions[s.id] !== undefined ? `: ${reactions[s.id]} ${outputUnit}` : ''}
                  </span>
                </p>
              </div>
            ))}
          </div>
        </div>

        {/* Diagram Gaya Geser (SFD) */}
        <div className="p-6 bg-gray-50 rounded-xl shadow-md border border-gray-200 mb-10">
          <h2 className="text-2xl font-semibold text-gray-700 mb-4 text-center">Diagram Gaya Geser (SFD)</h2>
          <ResponsiveContainer width="100%" height={300}>
            <LineChart
              data={sfdData}
              margin={{ top: 20, right: 30, left: 20, bottom: 5 }}
            >
              {/* SVG Definitions for Patterns */}
              <defs>
                {/* Pattern for positive shear (light blue vertical lines) */}
                <pattern id="pattern-positive-shear" x="0" y="0" width="10" height="10" patternUnits="userSpaceOnUse">
                  <line x1="0" y1="0" x2="0" y2="10" stroke="#AEC6CF" strokeWidth="1" /> {/* Lighter blue */}
                  <line x1="5" y1="0" x2="5" y2="10" stroke="#AEC6CF" strokeWidth="1" />
                </pattern>
                {/* Pattern for negative shear (darker blue vertical lines) */}
                <pattern id="pattern-negative-shear" x="0" y="0" width="10" height="10" patternUnits="userSpaceOnUse">
                  <line x1="0" y1="0" x2="0" y2="10" stroke="#6495ED" strokeWidth="1" /> {/* Cornflower blue, darker */}
                  <line x1="5" y1="0" x2="5" y2="10" stroke="#6495ED" strokeWidth="1" />
                </pattern>
              </defs>

              <CartesianGrid strokeDasharray="3 3" />
              <XAxis type="number" dataKey="x" domain={[0, beamLength]} label={{ value: "Posisi (m)", position: "insideBottomRight", offset: 0 }} ticks={generateXAxisTicks()} />
              <YAxis label={{ value: sfdYAxisLabel, angle: -90, position: "insideLeft" }} domain={sfdYAxisDomain} ticks={sfdYAxisTicks} />
              <Tooltip formatter={(value) => [`${value} ${outputUnit}`, 'Gaya Geser']} />
              <ReferenceLine y={0} stroke="#9ca3af" strokeDasharray="3 3" />

              {/* Area components for shading. V_pos and V_neg are created in sfdData */}
              <Area type="linear" dataKey="V_pos" stroke={false} fill="url(#pattern-positive-shear)" fillOpacity={0.8} connectNulls={false} />
              <Area type="linear" dataKey="V_neg" stroke={false} fill="url(#pattern-negative-shear)" fillOpacity={0.8} connectNulls={false} />

              {/* The actual SFD line */}
              <Line type="linear" dataKey="V" stroke="#4F46E5" strokeWidth={2} dot={false} />

              {/* LabelList for +/- signs */}
              <LabelList
                dataKey="V"
                position="outside" // Try 'outside' for better visibility
                content={({ x, y, value, index }) => {
                  // Only add label if value is not zero and at an "integer-like" position
                  // This prevents labels on every epsilon point
                  if (Math.abs(value) > 1e-9 && Math.abs(sfdData[index].x - Math.round(sfdData[index].x)) < epsilonForPlotting * 2) {
                    const sign = value > 0 ? '+' : '-';
                    const fill = value > 0 ? '#4F46E5' : '#CC0000'; // Blue for positive, Red for negative
                    const adjustedY = value > 0 ? y - 10 : y + 20; // Adjust y based on sign for better placement

                    return (
                      <text
                        x={x}
                        y={adjustedY}
                        fill={fill}
                        fontSize={14}
                        textAnchor="middle"
                        dominantBaseline="middle"
                        fontWeight="bold"
                      >
                        {sign}
                      </text>
                    );
                  }
                  return null;
                }}
              />
            </LineChart>
          </ResponsiveContainer>
        </div>

        {/* Diagram Momen Lentur (BMD) */}
        <div className="p-6 bg-gray-50 rounded-xl shadow-md border border-gray-200 mb-10">
          <h2 className="text-2xl font-semibold text-gray-700 mb-4 text-center">Diagram Momen Lentur (BMD)</h2>
          <ResponsiveContainer width="100%" height={300}>
            <LineChart
              data={bmdData} // Menggunakan bmdData (lebih detail)
              margin={{ top: 20, right: 30, left: 20, bottom: 5 }}
            >
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis type="number" dataKey="x" domain={[0, beamLength]} label={{ value: "Posisi (m)", position: "insideBottomRight", offset: 0 }} ticks={generateXAxisTicks()} />
              <YAxis label={{ value: bmdYAxisLabel, angle: -90, position: "insideLeft" }} domain={bmdYAxisDomain} ticks={bmdYAxisTicks} />
              <Tooltip formatter={(value) => [`${value} ${outputUnit === 'kN' ? 'kNm' : outputUnit + 'm'}`, 'Momen']}/>
              <Line type="linear" dataKey="M" stroke="#10B981" strokeWidth={2} dot={false} />
              <ReferenceLine y={0} stroke="#9ca3af" strokeDasharray="3 3" />
            </LineChart>
          </ResponsiveContainer>
        </div>

        {/* Tabel Hasil SFD */}
        <div className="p-6 bg-blue-50 rounded-xl shadow-md border border-blue-200 mb-10">
            <h2 className="text-2xl font-semibold text-blue-600 mb-4 text-center">Tabel Gaya Geser (SFD)</h2>
            <div className="overflow-x-auto">
                <table className="min-w-full bg-white rounded-lg shadow-sm text-sm">
                    <thead>
                        <tr className="bg-blue-100 text-blue-800 uppercase text-left">
                            <th className="py-2 px-4 border-b border-blue-200">Posisi (m)</th>
                            <th className="py-2 px-4 border-b border-blue-200">Gaya Geser ({outputUnit})</th>
                        </tr>
                    </thead>
                    <tbody>
                        {tableSfdData.map((data, index) => (
                            <tr key={index} className="hover:bg-blue-50">
                                <td className="py-2 px-4 border-b border-blue-200 font-medium">{data.x}</td>
                                <td className="py-2 px-4 border-b border-blue-200">{data.V}</td>
                            </tr>
                        ))}
                    </tbody>
                </table>
            </div>
        </div>

        {/* Tabel Hasil BMD */}
        <div className="p-6 bg-purple-50 rounded-xl shadow-md border border-purple-200">
            <h2 className="text-2xl font-semibold text-purple-600 mb-4 text-center">Tabel Momen Lentur (BMD)</h2>
            <div className="overflow-x-auto">
                <table className="min-w-full bg-white rounded-lg shadow-sm text-sm">
                    <thead>
                        <tr className="bg-purple-100 text-purple-800 uppercase text-left">
                            <th className="py-2 px-4 border-b border-purple-200">Posisi (m)</th>
                            <th className="py-2 px-4 border-b border-purple-200">Momen ({outputUnit === 'kN' ? 'kNm' : outputUnit + 'm'})</th>
                        </tr>
                    </thead>
                    <tbody>
                        {tableBmdData.map((data, index) => (
                            <tr key={index} className="hover:bg-purple-50">
                                <td className="py-2 px-4 border-b border-purple-200 font-medium">{data.x}</td>
                                <td className="py-2 px-4 border-b border-purple-200">{data.M}</td>
                            </tr>
                        ))}
                    </tbody>
                </table>
            </div>
        </div>

      </div>
    </div>
  );
}

export default App;
