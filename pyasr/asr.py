import os
import shutil
import subprocess
import pkg_resources

from Bio.Phylo.PAML import codeml

from .read import read_codeml_output

# paml 실행
def reconstruct(
    df,
    id_col='uid',
    sequence_col='sequence',
    working_dir='',
    save_ancestors=False,
    altall_cutoff=0.2,
    infer_gaps=True,
    aaRatefile='lg',
    **kwargs
    ):

    df = df.copy()

    # 기본 인수 구성
    default_options = dict(verbose=9, CodonFreq=None, cleandata=0,
        fix_blength=2, NSsites=None, fix_omega=None, clock=None,
        ncatG=8, runmode=0, fix_kappa=None, fix_alpha=1, Small_Diff=1.0e-6,
        method=0, Malpha=None, aaDist=None, RateAncestor=2, icode=None,
        alpha=None, seqtype=2, omega=None, getSE=None, noisy=3, Mgene=None,
        kappa=None, model=3, ndata=None)

    # 기본 인수를 현재 값으로 업데이트합니다.
    default_options.update(**kwargs)

    # ---------------- 모델 준비 ----------------
    # 패키지에서 프로젝트 디렉토리로 모델을 복사합니다.
    path_to_model = pkg_resources.resource_filename(
        'pyasr', os.path.join('dat', '{}.dat'.format(aaRatefile)))

    model_file = '{}.dat'.format(aaRatefile)
    model_path = os.path.join(working_dir, model_file)
    shutil.copyfile(path_to_model, model_path)

    # ----------------------

    curr_path = os.getcwd()
    proj_path = os.path.join(curr_path, working_dir)
    ali_path = os.path.join(working_dir, 'ali-to-reconstruct.phy')
    tree_path = os.path.join(working_dir, 'tree-to-reconstruct.phy')
    out_path = os.path.join(working_dir, 'results.txt')
    ctl_path = os.path.join(working_dir, 'codeml_options.ctl')
    rst_path = os.path.join(working_dir, 'rst')

    df.phylo.to_fasta(
        filename=ali_path,
        id_col=id_col,
        sequence_col=sequence_col,
    )

    df.phylo.to_newick(
        filename=tree_path,
        taxon_col=id_col,
        node_col=id_col,
        suppress_internal_node_labels=True,
    )

    df.phylo.to_newick(
        taxon_col=id_col,
        node_col=id_col,
        suppress_internal_node_labels=True,
    )

    # 제어 파일을 빌드하고 저장합니다.
    cml = codeml.Codeml(alignment=ali_path,
        tree=tree_path,
        out_file=out_path,
        working_dir=working_dir)
    cml.set_options(aaRatefile=model_file, **default_options)
    cml.ctl_file = ctl_path
    cml.write_ctl_file()

    # ----------------------

    os.chdir(proj_path)
    output = subprocess.run(['codeml', 'codeml_options.ctl'])
    os.chdir(curr_path)

    # ----------------------

    return read_codeml_output(rst_path, df)
