import dendropy

def infer_gaps_in_tree(df_seq, tree, id_col='id', sequence_col='sequence'):
    """문자 매트릭스를 DendroPy 트리에 추가하고 Fitch 알고리즘을 사용하여 갭을 추론합니다.

    조상 노드의 서열에서 갭을 추론합니다.
    """
    taxa = tree.taxon_namespace

    # fasta로 정렬을 가져옵니다.
    alignment = df_seq.phylo.to_fasta(id_col=id_col, id_only=True,
                                      sequence_col=sequence_col)

    # Dendropy에서 서열 데이터 매트릭스를 빌드합니다.
    data = dendropy.ProteinCharacterMatrix.get(
        data=alignment,
        schema="fasta",
        taxon_namespace=taxa)

    # 서열 데이터와 트리 데이터 간의 맵 객체를 구성합니다.
    taxon_state_sets_map = data.taxon_state_sets_map(gaps_as_missing=False)

    # 갭 배치를 결정하는 Fitch 알고리즘
    dendropy.model.parsimony.fitch_down_pass(tree.postorder_node_iter(),
            taxon_state_sets_map=taxon_state_sets_map)
    dendropy.model.parsimony.fitch_up_pass(tree.preorder_node_iter())
    return tree
