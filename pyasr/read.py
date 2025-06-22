import re
import dendropy
from phylopandas import DataFrame


def read_codeml_output(
    filename,
    df,
    altall_cutoff=0.2,
    ):
    """codeml 파일을 읽습니다.
    """
    # paml 출력을 읽습니다.
    with open(filename, 'r') as f:
        data = f.read()

    # codeml 출력에서 모든 트리를 추출합니다.
    regex = re.compile('\([()\w\:. ,]+;')
    trees = regex.findall(data)
    anc_tree = trees[2]

    # codeml 파일의 첫 번째 트리는 원본 입력 트리입니다.
    tip_tree = dendropy.Tree.get(data=trees[0], schema='newick')

    # codeml 파일의 세 번째 트리는 조상 트리입니다.
    anc_tree = dendropy.Tree.get(data=trees[2], schema='newick')

    # 반환할 메인 트리
    tree = tip_tree

    # 조상을 메인 트리 객체에 매핑합니다.
    ancestors = anc_tree.internal_nodes()
    for i, node in enumerate(tree.internal_nodes()):
        node.label = ancestors[i].label

    # 노드를 데이터프레임에 매핑합니다.
    df['reconstruct_label'] = None
    for node in tree.postorder_node_iter():

        # 부모 노드를 무시합니다.
        if node.parent_node is None:
            pass

        elif node.is_leaf():
            node_label = node.taxon.label
            parent_label = node.parent_node.label
            # 노드 레이블을 설정합니다.
            df.loc[df.uid == node_label, 'reconstruct_label'] = node_label

            # 부모 레이블을 설정합니다.
            parent_id = df.loc[df.uid == node_label, 'parent'].values[0]
            df.loc[df.id == parent_id, 'reconstruct_label'] = node.parent_node.label

        elif node.is_internal():
            label = node.label
            parent_id = df.loc[df.reconstruct_label == label, 'parent'].values[0]
            df.loc[df.id == parent_id, 'reconstruct_label'] = node.parent_node.label


    # 내부 노드에 대한 데이터 블록을 찾기 위한 정규 표현식을 컴파일합니다.
    node_regex = re.compile("""Prob distribution at node [0-9]+, by site[-\w():.\s]+\n""")

    # 이 데이터 블록에서 노드 번호를 추출합니다.
    node_num_regex = re.compile("[0-9]+")

    # 모든 조상에 대한 데이터프레임을 가져옵니다.
    df['ml_sequence'] = None
    df['ml_posterior'] = None
    df['alt_sequence'] = None
    df['alt_posterior'] = None

    for node in node_regex.findall(data):
        # 노드 레이블 가져오기
        node_label = node_num_regex.search(node).group(0)

        # 사이트 데이터와 일치하는 정규식 컴파일
        site_regex = re.compile("(?:\w\(\w.\w{3}\) )+")

        # 사이트 데이터의 각 일치 항목을 반복합니다.
        ml_sequence, ml_posterior, alt_sequence, alt_posterior = [], [], [], []

        for site in site_regex.findall(node):

            # 잔기를 반복합니다.
            scores = [float(site[i+2:i+7]) for i in range(0,len(site), 9)]
            residues = [site[i] for i in range(0, len(site), 9)]

            # 정렬된 점수의 인덱스를 가져옵니다.
            sorted_score_index = [i[0] for i in sorted(
                enumerate(scores),
                key=lambda x:x[1],
                reverse=True)]

            ml_idx = sorted_score_index[0]
            alt_idx = sorted_score_index[1]

            # 대체 사이트를 유지해야 하는지 확인합니다.
            ml_sequence.append(residues[ml_idx])
            ml_posterior.append(scores[ml_idx])

            if scores[alt_idx] < altall_cutoff:
                alt_idx = ml_idx

            alt_sequence.append(residues[alt_idx])
            alt_posterior.append(scores[alt_idx])

        keys = [
            "ml_sequence",
            "ml_posterior",
            "alt_sequence",
            "alt_posterior"
        ]

        vals = [
            "".join(ml_sequence),
            sum(ml_posterior) / len(ml_posterior),
            "".join(alt_sequence),
            sum(alt_posterior) / len(alt_posterior),
        ]

        df.loc[df.reconstruct_label == node_label, keys] = vals

    return df
