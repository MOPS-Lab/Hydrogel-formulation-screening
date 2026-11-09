/* Input uncertainty only: conditional on the fixed GP and on gel formation.
Uniform scenario ranges are user assumptions, not measured laboratory errors. */
function createMonteCarlo(model,meta){
 function rng(seed){let x=seed>>>0;return ()=>{x=(x+0x6D2B79F5)>>>0;let t=x;t=Math.imul(t^(t>>>15),t|1);t^=t+Math.imul(t^(t>>>7),t|61);return ((t^(t>>>14))>>>0)/4294967296;};}
 const quantile=(x,p)=>{const k=(x.length-1)*p,i=Math.floor(k);return x[i]+(x[Math.min(i+1,x.length-1)]-x[i])*(k-i);};
 function config(o={}){
  const q={n:o.n??5000,seed:o.seed??20261014,stock_halfwidth_pp:o.stock_halfwidth_pp??.1,volume_halfwidth_percent:o.volume_halfwidth_percent??2,crosslinker_halfwidth_percent:o.crosslinker_halfwidth_percent??5,target_Pa:o.target_Pa??15000,tolerance_percent:o.tolerance_percent??10};
  if(!Number.isInteger(q.n)||q.n<100||q.n>50000||!Number.isInteger(q.seed)||q.seed<0||q.seed>4294967295)throw Error('N: celo število 100–50000; seme: celo število 0–4294967295.');
  if(!Number.isFinite(q.stock_halfwidth_pp)||q.stock_halfwidth_pp<0||q.stock_halfwidth_pp>1||!Number.isFinite(q.volume_halfwidth_percent)||q.volume_halfwidth_percent<0||q.volume_halfwidth_percent>=100||!Number.isFinite(q.crosslinker_halfwidth_percent)||q.crosslinker_halfwidth_percent<0||q.crosslinker_halfwidth_percent>=100)throw Error('Neveljavni razponi: koncentracija 0–1 odstotne točke, relativna razpona 0–manj kot 100 %.');
  if(!Number.isFinite(q.target_Pa)||q.target_Pa<=0||!Number.isFinite(q.tolerance_percent)||q.tolerance_percent<=0||q.tolerance_percent>100)throw Error('Ciljni G mora biti pozitiven; toleranca 0–100 %.');return q;
 }
 function drawRecipe(o,q,random){
  const u=Array.from({length:8},()=>2*random()-1),stocks=o.stocks.map((v,i)=>v+q.stock_halfwidth_pp*u[i]),vol=o.fractions.map((v,i)=>v*(1+q.volume_halfwidth_percent*u[3+i]/100)),vx=1+q.volume_halfwidth_percent*u[6]/100,total=vol.reduce((a,b)=>a+b,0)+vx,ionstock=o.dose*(1+q.crosslinker_halfwidth_percent*u[7]/100);
  if(stocks.some((s,i)=>o.fractions[i]>0&&s<=0)||total<=0||ionstock<0)return null;
  return {c:stocks.map((v,i)=>v*vol[i]/total),dose:ionstock*vx/total,ratio:vol.reduce((a,b)=>a+b,0)/vx};
 }
 async function simulate(o,options={},progress=()=>{}){
  const q=config(options),nominal=model.predict(o);if(nominal.G_Pa===null)return {status:'observed_no_gel',nominal,config:q};
  const ion=nominal.ion,random=rng(q.seed),values=[],ratios=[],outside=[],failed=meta.records.filter(r=>!r.gel_formed&&r.ion===ion);let invalid=0,noGel=0,below=0,above=0,hits=0;const doseRange=meta.support[ion];
  for(let i=0;i<q.n;i++){
   const r=drawRecipe(o,q,random);if(!r){invalid++;continue;}
   if(failed.some(f=>Math.abs(f.dose_M-r.dose)<1e-8&&f.c.every((v,j)=>Math.abs(v-r.c[j])<1e-8))){noGel++;continue;}
   const G=10**model.mean(r.c,ion,r.dose);if(!Number.isFinite(G)||G<=0){invalid++;continue;}
   values.push(G);ratios.push(r.ratio);outside.push(!model.inside(r.c,ion,r.dose));if(r.dose<doseRange.dose_min)below++;if(r.dose>doseRange.dose_max)above++;if(Math.abs(G/q.target_Pa-1)<=q.tolerance_percent/100)hits++;
   if((i+1)%250===0){progress({done:i+1,total:q.n});await new Promise(resolve=>setTimeout(resolve,0));}
  }
  if(!values.length)return {status:'no_valid_simulations',nominal,config:q,invalid,noGel};
  values.sort((a,b)=>a-b);ratios.sort((a,b)=>a-b);const mean=values.reduce((a,b)=>a+b,0)/values.length,sd=Math.sqrt(values.reduce((s,v)=>s+(v-mean)**2,0)/Math.max(1,values.length-1)),p=hits/values.length,lo=values[0],hi=values.at(-1),bins=Array(20).fill(0);
  values.forEach(v=>bins[Math.min(19,Math.floor((v-lo)/Math.max(hi-lo,1e-12)*20))]++);
  const z=1.959963984540054,den=1+z*z/values.length,center=(p+z*z/(2*values.length))/den,half=z*Math.sqrt(p*(1-p)/values.length+z*z/(4*values.length**2))/den;
  return {status:'ok',nominal,config:q,valid_n:values.length,invalid_n:invalid,known_no_gel_n:noGel,mean_Pa:mean,median_Pa:quantile(values,.5),sd_Pa:sd,preparation90_Pa:[quantile(values,.05),quantile(values,.95)],target_fraction_given_valid:p,target_fraction_MC_se:Math.sqrt(p*(1-p)/values.length),target_fraction_numerical95:[Math.max(0,center-half),Math.min(1,center+half)],p90_relative_target_error:quantile(values.map(v=>Math.abs(v/q.target_Pa-1)).sort((a,b)=>a-b),.9),outside_training_fraction:outside.filter(Boolean).length/values.length,outside_dose_count:below+above,polymer_to_crosslinker_volume_ratio90:[quantile(ratios,.05),quantile(ratios,.95)],histogram:{min_Pa:lo,max_Pa:hi,counts:bins},conditional_on_gel_formation:true,includes_model_uncertainty:false,distribution:'independent uniform scenarios on stock concentrations and dispensed volumes',values_Pa:options.include_values?values:undefined};
 }
 async function robustInverse(o={},options={},progress=()=>{}){
  const q=config({...options,target_Pa:o.target_Pa,tolerance_percent:o.tolerance_percent??5});
  const first=await model.inverse({...o,max_recipes:24},p=>progress({stage:'grid',...p}));let list=[];
  for(let i=0;i<first.recipes.length;i++){const r=first.recipes[i],p=r.prediction,recipe={stocks:p.stocks_percent,fractions:p.liquid_fractions,ion:p.ion,dose:p.crosslinker_stock_M};const mc=await simulate(recipe,{...q,n:Math.min(1000,q.n)});if(mc.status==='ok')list.push({...r,mc,recipe});progress({stage:'screen',done:i+1,total:first.recipes.length});}
  const order=(a,b)=>Number(b.within_tolerance)-Number(a.within_tolerance)||Number(b.prediction.inside_combined_input_domain)-Number(a.prediction.inside_combined_input_domain)||b.mc.target_fraction_given_valid-a.mc.target_fraction_given_valid||a.mc.p90_relative_target_error-b.mc.p90_relative_target_error;
  list.sort(order);list=list.slice(0,8);
  for(let i=0;i<list.length;i++){list[i].mc=await simulate(list[i].recipe,q);progress({stage:'refine',done:i+1,total:list.length});}list.sort(order);
  return {target_Pa:o.target_Pa,recipes:list,examined:first.examined,MC_shortlist_count:first.recipes.length,config:q,global_optimum_guaranteed:false,ranking:'nominal tolerance, nominal measured domain, scenario target fraction, 90th-percentile target error; common random numbers'};
 }
 return {simulate,robustInverse,drawRecipe,rng,config};
}
