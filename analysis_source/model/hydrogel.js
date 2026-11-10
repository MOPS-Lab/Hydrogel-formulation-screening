/* v4: source-aware GP; no measured target is a prediction-time input. */
function createHydrogelModel(m,meta,diag,alternative){
 const names=['ALG','TOCNF','SCLG'],family=c=>names.filter((_,i)=>c[i]>1e-10).join('+');
 const failed=meta.records.filter(r=>!r.gel_formed),positive=meta.records.filter(r=>r.gel_formed);
 function sourceFor(ion){return ion==='none'?'legacy':'KS';}
 function mean(c,ion,dose,source=sourceFor(ion)){
  const x=[...c,Math.log10(1+dose/.005)].map((v,i)=>(v-m.scaler_mean[i])/m.scaler_scale[i]);let s=0;
  for(let j=0;j<m.X.length;j++){let d=0;for(let k=0;k<4;k++)d+=(x[k]-m.X[j][k])**2;const t=Math.sqrt(5*d)/m.length;s+=(1+t+t*t/3)*Math.exp(-t)*(ion===m.ions[j]?1:m.rho_ion)*(source===m.sources[j]?1:m.rho_source)*m.alpha[j];}
  return m.y_mean+s;
 }
 function alternativeLog(c,ion,dose,source=sourceFor(ion)){
  const C=c.reduce((a,b)=>a+b,0),w=c.map(v=>v/C),d=Math.log10(1+dose/.005),lc=Math.log10(C),is=['Ca','Mg','Ba','Sr','Fe','Cu','Al','none'].map(x=>Number(x===ion)),src=Number(source==='KS');
  const x=[lc,d,w[0],w[1],...is,src,...is.map(v=>v*src),lc*w[0],lc*w[1],d*w[0],d*w[1]];
  return alternative.intercept+x.reduce((a,v,i)=>a+(v-alternative.scaler_mean[i])/alternative.scaler_scale[i]*alternative.coef[i],0);
 }
 function evidence(c,ion,dose,G){
  const fam=family(c),count=positive.filter(r=>r.source==='KS'&&r.ion===ion&&family(r.c)===fam).length;
  let v=diag.by_family_ion[fam+'|'+ion],scope='ista družina in ion';
  if(!v||v.n<5||v.n_compositions<3){v=diag.by_family[fam];scope='ista družina, različni ioni';}
  if(!v||v.n<5||v.n_compositions<3){v=diag.all;scope='celoten nabor novih sestav';}
  const supportedStats=count>0?{...v,scope,reference_range_Pa:[G/v.factor80,G*v.factor80],definition:'Descriptive 80th percentile of held-out errors, not a calibrated prediction interval'}:null;
  const alt=10**alternativeLog(c,ion,dose),disagreement=Math.max(G/alt,alt/G);
  return {validation_reference:supportedStats,alternative_G_Pa:alt,model_disagreement_factor:disagreement,positive_KS_family_ion_count:count};
 }
 function sensitivity(o){
  const p=predict(o);if(p.G_Pa===null)return null;
  const stocks=o.stocks.slice(),fractions=o.fractions.slice(),dose=o.dose,ion=o.ion;
  const cases=[];const add=(ss,ff,dd)=>{try{const cc=ss.map((v,i)=>v*ff[i]/2),q=meta.support[ion];if(dd/2<q.dose_min||dd/2>q.dose_max||failures(cc,ion,dd/2).length)return;cases.push({G:10**mean(cc,ion,dd/2),inside:inside(cc,ion,dd/2)});}catch{}};
  for(let i=0;i<3;i++)if(fractions[i]>0)for(const sign of [-1,1]){const ss=stocks.slice();ss[i]+=.1*sign;if(ss[i]>=.1&&ss[i]<=3)add(ss,fractions,dose);}
  for(let i=0;i<3;i++)for(let j=i+1;j<3;j++)if(fractions[i]>0&&fractions[j]>0)for(const sign of [-1,1]){const ff=fractions.slice();ff[i]+=.02*sign;ff[j]-=.02*sign;if(ff.every(v=>v>=0&&v<=1))add(stocks,ff,dose);}
  for(const k of [.95,1.05])add(stocks,fractions,dose*k);
  if(!cases.length)return null;const vals=[p.G_Pa,...cases.map(x=>x.G)];return {min_Pa:Math.min(...vals),max_Pa:Math.max(...vals),max_change_percent:100*Math.max(...vals.map(v=>Math.abs(v/p.G_Pa-1))),scenarios:cases.length,outside_domain_count:cases.filter(x=>!x.inside).length,definition:'One input at a time: stock +/-0.1 percentage point; pair fractions +/-0.02; crosslinker stock +/-5%; model scenario, not a predictive interval'};
 }
 function inside(c,ion,dose,domains=meta.domains){
  const q=domains[family(c)+'|'+ion];if(!q)return false;
  const d=[...c,Math.log10(1+dose/.005)].map((v,i)=>v-q.origin[i]);if(!q.rank)return Math.hypot(...d)<1e-7;
  const z=Array.from({length:q.rank},(_,j)=>d.reduce((s,v,i)=>s+v*q.basis[i][j],0));
  if(Math.hypot(...d.map((v,i)=>v-q.basis[i].reduce((s,u,j)=>s+u*z[j],0)))>1e-7)return false;
  if(q.rank===1)return z[0]>=q.lower-1e-7&&z[0]<=q.upper+1e-7;
  return q.equations.every(e=>e.slice(0,-1).reduce((s,u,j)=>s+u*z[j],e.at(-1))<=1e-7);
 }
 function failures(c,ion,dose){return failed.filter(r=>r.ion===ion&&Math.abs(r.dose_M-dose)<1e-8&&r.c.every((v,i)=>Math.abs(v-c[i])<1e-8)).map(r=>r.id);}
 function predictEffective(c,ion='Ca',dose=.045){
  if(!Array.isArray(c)||c.length!==3||!c.every(v=>Number.isFinite(v)&&v>=0)||c.reduce((s,v)=>s+v,0)<=0||!Number.isFinite(dose)||dose<0)throw Error('Neveljavna sestava ali odmerek.');
  if(!Object.hasOwn(meta.support,ion))throw Error('Neznan ion.');if(ion==='none'&&dose!==0)throw Error('Brez ionov zahteva odmerek 0.');if(dose===0)ion='none';
  const C=c.reduce((s,v)=>s+v,0),fam=family(c),q=meta.support[ion];
  if(C>1.5+1e-8)throw Error('Skupna koncentracija nad 1,5 % presega obseg kalkulatorja.');
  if(dose<q.dose_min-1e-9||dose>q.dose_max+1e-9)throw Error('Odmerek zunaj razpona podatkov za '+ion+'.');
  const inData=inside(c,ion,dose),inKS=inside(c,ion,dose,meta.ks_domains),ids=failures(c,ion,dose);
  const common={version:meta.version,c_percent:c.slice(),total_percent:C,polymer_mass_fractions:c.map(v=>v/C),family:fam,ion,dose_M:dose,prediction_source:sourceFor(ion),inside_combined_input_domain:inData,inside_KS_input_domain:inKS,gel_probability:null,prediction_interval:null,conditional_on_gel_formation:true};
  if(ids.length)return {...common,status:'observed_no_gel',G_Pa:null,known_failure_ids:ids,flags:['observed_failure_not_classifier_prediction']};
  if(ion==='none'&&!inData)throw Error('Brez ionov je dovoljena le napoved v območju meritev.');
  const flags=['exploratory_model'];if(!inData)flags.push('input_extrapolation');if(!inKS)flags.push('outside_KS_input_domain');if(C<.5)flags.push('total_below_inverse_range');
  const count=positive.filter(r=>r.source==='KS'&&r.ion===ion&&family(r.c)===fam).length;
  if(!count)flags.push('no_KS_family_ion_measurements');if(ion==='Mg')flags.push('no_KS_Mg_validation');
  const relevant=positive.filter(r=>r.ion===ion).map(r=>({id:r.id,source:r.source,c:r.c,dose_M:r.dose_M,G_Pa:r.G_observed_Pa,distance:Math.hypot(...[...r.c,Math.log10(1+r.dose_M/.005)].map((v,i)=>(v-[...c,Math.log10(1+dose/.005)][i])/m.scaler_scale[i]))})).sort((a,b)=>a.distance-b.distance).slice(0,3);
  const G=10**mean(c,ion,dose);return {...common,status:'conditional_estimate',G_Pa:G,KS_family_ion_positive_count:count,nearest_observations:relevant,flags,...evidence(c,ion,dose,G)};
 }
 function predict(o={}){
  const stocks=o.stocks??[2,2,2],fractions=o.fractions??[.5,.5,0];
  if(stocks.length!==3||fractions.length!==3||!stocks.every(v=>Number.isFinite(v)&&v>=.1&&v<=3)||!fractions.every(v=>Number.isFinite(v)&&v>=0&&v<=1)||Math.abs(fractions.reduce((a,b)=>a+b,0)-1)>1e-8)throw Error('Zaloge: 0,1–3,0 %. Deleži: 0–1, vsota mora biti 1.');
  return {...predictEffective(stocks.map((v,i)=>v*fractions[i]/2),o.ion??'Ca',(o.dose??.09)/2),stocks_percent:stocks.slice(),liquid_fractions:fractions.slice(),crosslinker_stock_M:o.dose??.09};
 }
 async function inverse(o={},progress=()=>{}){
  const target=o.target_Pa,tol=o.tolerance_percent??5,wanted=o.family??'ALG+TOCNF',selected=o.ion??'Ca';
  if(!Number.isFinite(target)||target<=0||!Number.isFinite(tol)||tol<=0||tol>50)throw Error('Ciljni G mora biti pozitiven; toleranca večja od 0 in največ 50 %.');
  const fams=['ALG','TOCNF','SCLG','ALG+TOCNF','ALG+SCLG','TOCNF+SCLG','ALG+TOCNF+SCLG'];if(wanted!=='any'&&!fams.includes(wanted))throw Error('Neznana družina.');
  const ions=selected==='all'?Object.keys(meta.support).filter(x=>x!=='none'):[selected];if(ions.some(x=>!Object.hasOwn(meta.support,x)))throw Error('Neznan ion.');
  const ws=[];for(let a=0;a<=20;a++)for(let t=0;t<=20-a;t++){const w=[a/20,t/20,(20-a-t)/20];if(wanted==='any'||family(w)===wanted)ws.push(w);}
  const logT=Math.log10(target),seeds=[];let examined=0,done=0;
  for(const ion of ions){
   const q=meta.support[ion],ds=ion==='none'?[0]:[...new Set([...Array.from({length:15},(_,i)=>q.dose_min*(q.dose_max/q.dose_min)**(i/14)),...meta.records.filter(r=>r.source==='KS'&&r.ion===ion).map(r=>r.dose_M)])].sort((a,b)=>a-b);
   for(let ci=10;ci<=30;ci++){
    for(const w of ws){const C=ci/10,c=w.map(v=>v*C/2);let best=null;
     for(let j=0;j<ds.length;j++){const dose=ds[j];examined++;if(failures(c,ion,dose).length)continue;
      const supported=inside(c,ion,dose);if((o.supported_only||ion==='none')&&!supported)continue;
      const logG=mean(c,ion,dose),err=Math.abs(logG-logT);if(!best||err<best.err)best={C,w,c,ion,dose,j,ds,logG,err,supported};
     }if(best)seeds.push(best);
    }done++;progress({done,total:ions.length*21,examined});await new Promise(resolve=>setTimeout(resolve,0));
   }
  }
  const shortlist=[...seeds.filter(x=>x.supported).sort((a,b)=>a.err-b.err).slice(0,80),...seeds.filter(x=>!x.supported).sort((a,b)=>a.err-b.err).slice(0,120)],results=[];
  for(const seed of shortlist){let best={...seed};
   for(const k of [seed.j-1,seed.j+1]){if(seed.ion==='none'||k<0||k>=seed.ds.length)continue;
    let lo=Math.min(seed.dose,seed.ds[k]),hi=Math.max(seed.dose,seed.ds[k]),vl=mean(seed.c,seed.ion,lo)-logT,vh=mean(seed.c,seed.ion,hi)-logT;if(vl*vh>0)continue;
    for(let j=0;j<20;j++){const d=Math.sqrt(lo*hi),lg=mean(seed.c,seed.ion,d),v=lg-logT,sup=inside(seed.c,seed.ion,d);if(Math.abs(v)<best.err&&(!o.supported_only||sup)&&(!seed.supported||sup)&&!failures(seed.c,seed.ion,d).length)best={...seed,dose:d,logG:lg,err:Math.abs(v),supported:sup};if(vl*v<=0){hi=d;vh=v;}else{lo=d;vl=v;}}
   }
   const prediction=predict({stocks:[best.C,best.C,best.C],fractions:best.w,ion:best.ion,dose:2*best.dose});if(prediction.G_Pa===null)continue;
   const error_percent=100*(prediction.G_Pa/target-1);results.push({prediction,error_percent,within_tolerance:Math.abs(error_percent)<=tol,sensitivity:sensitivity({stocks:prediction.stocks_percent,fractions:prediction.liquid_fractions,ion:prediction.ion,dose:prediction.crosslinker_stock_M})});
  }
  results.sort((a,b)=>Number(b.within_tolerance)-Number(a.within_tolerance)||(a.within_tolerance&&b.within_tolerance?Number(b.prediction.inside_combined_input_domain)-Number(a.prediction.inside_combined_input_domain):0)||(a.within_tolerance&&b.within_tolerance?(a.sensitivity?.max_change_percent??Infinity)-(b.sensitivity?.max_change_percent??Infinity):0)||Math.abs(a.error_percent)-Math.abs(b.error_percent));
  const recipes=[];for(const r of results){const p=r.prediction;if(recipes.some(s=>s.prediction.ion===p.ion&&Math.abs(s.prediction.total_percent-p.total_percent)<.1495&&s.prediction.polymer_mass_fractions.reduce((v,w,i)=>v+Math.abs(w-p.polymer_mass_fractions[i]),0)<.149))continue;recipes.push(r);if(recipes.length===(o.max_recipes??8))break;}
  return {target_Pa:target,tolerance_percent:tol,recipes,examined,search:{stocks_percent:[1,3],stock_step:.1,fraction_step:.05,dose:'15 logarithmic points plus KS doses, local refinement',canonical_stocks:'equal stock concentrations for active polymers',global_optimum_guaranteed:false},conditional_on_gel_formation:true,prediction_interval:null};
 }
 return {predict,predictEffective,inverse,mean,inside,family,sensitivity,alternativeLog};
}
