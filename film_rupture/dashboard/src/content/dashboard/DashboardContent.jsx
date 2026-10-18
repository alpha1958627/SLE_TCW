import React from "react";
import { DataComponent, DataTable, EvidenceChart, Filters, SectionHeader,
  SortableItem, SortableRegion, useDataApp } from "../../data-app-public.jsx";
import "./simulation.css";

const palette = {"run 1":"#202020","run 2":"#be3434","run 3":"#2a62b3","run 4":"#26816f","run 5":"#ad750a","old r1":"#7b3fa2"};
const f = (value, digits=3) => Number.isFinite(value) ? value.toFixed(digits) : "—";
const percent = value => Number.isFinite(value) ? `${value > 0 ? "+" : ""}${value.toFixed(3)}%` : "—";
const columns = [
  {field:"ratio",label:"厚度倍率"}, {field:"curve",label:"run / 对照"},
  {field:"server",label:"服务器"}, {field:"state",label:"状态"},
  {field:"step",label:"step"}, {field:"t",label:"t",renderCell:v=>f(v,1)},
  {field:"hmin",label:"hmin",renderCell:v=>f(v,4)},
  {field:"logDriftPct",label:"V 最新漂移",renderCell:v=>percent(v)},
  {field:"maxLogDriftPct",label:"V 最大偏差",renderCell:v=>percent(v)},
  {field:"fieldDriftPct",label:"最近全场 V 偏差",renderCell:v=>percent(v)},
  {field:"job",label:"作业"},
];

export function DashboardContent() {
  const { snapshot, queries, filters, setFilter, reviewedRows } = useDataApp();
  const runs = reviewedRows("run_status", ["task"]);
  const curves = reviewedRows("film_history", ["task","t","curve"]);
  const volumes = reviewedRows("field_volume", ["task","t"]);
  const maxDrift = runs.length ? Math.max(...runs.map(r=>r.maxLogDriftPct)) : null;
  const active = runs.filter(r=>r.state === "RUNNING").length;
  const alerts = runs.filter(r=>r.maxLogDriftPct > 0.1).length;
  const groups = ["1.05","1.10","1.15","1.20"];
  const rows = [
    {id:"simulation:smallmultiples1",items:["hmin-h105","hmin-h110"]},
    {id:"simulation:smallmultiples2",items:["hmin-h115","hmin-h120"]},
    {id:"simulation:volume",items:["volume-check"]},
    {id:"simulation:records",items:["simulation-records"]},
  ];
  return <article className="simulation-dashboard">
    <Filters sticky filters={snapshot.filters ?? []} queries={queries} values={filters} onChange={setFilter}
      trailingControls={<span className="sim-cadence">每 12 小时更新 · 北京时间</span>} />
    <DataComponent id="simulation-summary" queryId="run_status" kind="metrics" title="本批轨迹" displayRows={runs} sourceRows={runs} variant="plain">
      <div className="sim-summary" data-reviewed-rows>
        <span><strong>{runs.length}</strong> 条已选轨迹</span>
        <span><strong>{active}</strong> 运行中</span>
        <span><strong>{alerts}</strong> 条 V 偏差超过 0.1%</span>
        <span>最大日志相对偏差 <strong>{f(maxDrift,3)}%</strong></span>
      </div>
      <p className="sim-note">数据快照：{new Date(snapshot.generatedAt).toLocaleString("zh-CN",{timeZone:"Asia/Shanghai",hour12:false})}（北京时间）；每 12 小时更新。</p><p className="sim-note">t 为各工况的原始数值时间；h 已按初始膜厚归一化。曲线止于各 run 最新记录，尚未结束的轨迹不计作破碎寿命。</p>
    </DataComponent>
    <SectionHeader id="hmin-curves-heading" title="hmin vs t" />
    <SortableRegion id="simulation:canvas" label="液膜曲线和体积监测" variant="canvas" spacing="standard" columns={12} authoredRevision={1} rows={rows}>
      {groups.map(ratio=>{
        const id = "hmin-h" + ratio.replace(".","");
        const selected = curves.filter(r=>r.ratio===ratio);
        return <SortableItem key={id} id={id} label={`${ratio} h0`} kind="chart" span={6}>
          <EvidenceChart id={id} queryId="film_history" title={`${ratio} h₀ · hmin–t`}
            spec={{type:"line",x:"t",y:"hmin",series:"curve",colors:palette,dashedFields:["old r1"],
              stackable:false,xLabel:"t（数值时间）",yLabel:"hmin / h初始",valueDecimals:4,yTickCount:6}}
            rows={selected} sourceRows={selected} variant="card" height={300}>
            <p className="sim-note">{selected.length ? `${new Set(selected.map(r=>r.task)).size} 条轨迹` : "当前筛选无轨迹"}{ratio==="1.05" ? "；紫色虚线 old r1 与 run 1 同种子。" : "；run 1–5 使用五个固定种子。"}</p>
          </EvidenceChart>
        </SortableItem>;
      })}
      <SortableItem id="volume-check" label="体积漂移" kind="chart" span={12}>
        <EvidenceChart id="volume-check" queryId="field_volume" title="独立全场积分：V 相对初始体积的偏差" variant="card"
          spec={{type:"line",x:"t",y:"driftPct",series:"identity",stackable:false,xLabel:"已保存全场的 t",yLabel:"偏差（%）",valueDecimals:3,startAtZero:false}}
          rows={volumes.map(r=>({...r,identity:`${r.ratio} · ${r.curve}`}))} sourceRows={volumes} height={280}>
          <p className="sim-alert" role="status">V 守恒尚未通过。日志原始 V 存在约 8 倍 MPI 重复归约；此图直接对保存的 64×64 全场积分，没有归一化场。全场时刻可能落后于最新日志。</p>
        </EvidenceChart>
      </SortableItem>
      <SortableItem id="simulation-records" label="逐 run 状态" kind="table" span={12}>
        <DataComponent id="simulation-records" queryId="run_status" title="逐 run 状态与数值检查" kind="table" variant="card" displayRows={runs} sourceRows={runs}>
          <DataTable rows={runs} columns={columns} rowKey="task" caption="21 条保留轨迹的状态、hmin 和体积漂移" />
          <p className="sim-note">V 最新 / 最大偏差来自日志相对首条值；全场偏差由独立积分得到，两者时刻可不同。0.1% 是监测警报线，不是已验证的物理验收标准。</p>
        </DataComponent>
      </SortableItem>
    </SortableRegion>
  </article>;
}