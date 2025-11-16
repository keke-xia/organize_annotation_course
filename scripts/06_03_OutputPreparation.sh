# 回到 MAKER 运行目录
cd /data/users/kxia/organize_annotation_course/gene_annotation

# 设定工具路径（课程自带的 MAKER 工具）
MAKERBIN=/data/courses/assembly-annotation-course/CDS_annotation/softwares/Maker_v3.01.03/src/bin

# 指定 master index 的绝对路径
IDX=/data/users/kxia/organize_annotation_course/gene_annotation/assembly.maker.output/assembly_master_datastore_index.log

# 合并 GFF（含序列）
$MAKERBIN/gff3_merge -s -d "$IDX" > assembly.all.maker.gff

# 合并 GFF（不含序列，后续处理常用这个）
$MAKERBIN/gff3_merge -n -s -d "$IDX" > assembly.all.maker.noseq.gff

# 合并蛋白/转录本 FASTA（会生成 assembly.all.maker.proteins.fasta / transcripts.fasta）
$MAKERBIN/fasta_merge -d "$IDX" -o assembly

# 看看产物
ls -lh assembly.all.maker.gff assembly.all.maker.noseq.gff assembly.all.maker.proteins.fasta assembly.all.maker.transcripts.fasta


cd "$WORKDIR"

# --- Inputs produced in Step 5 ---
protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"
gff="assembly.all.maker.noseq.gff"

# --- Prepare final dir and copies ---
mkdir -p final
cp -f "$gff"        "final/${gff}.renamed.gff"
cp -f "$protein"    "final/${protein}.renamed.fasta"
cp -f "$transcript" "final/${transcript}.renamed.fasta"

cd final

# --- Map IDs using MAKER utilities ---
# maker_map_ids creates an id.map; map_* apply it to GFF/FASTA
PREFIX="Pa-1"
"$MAKERBIN/maker_map_ids" --prefix "$PREFIX" --justify 7 "${gff}.renamed.gff" > id.map
"$MAKERBIN/map_gff_ids"   id.map "${gff}.renamed.gff"
"$MAKERBIN/map_fasta_ids" id.map "${protein}.renamed.fasta"
"$MAKERBIN/map_fasta_ids" id.map "${transcript}.renamed.fasta"

echo "[OK] Renamed with prefix: $PREFIX"