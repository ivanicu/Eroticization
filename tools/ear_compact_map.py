"""tools/ear_compact_map.py — the EAR arc merge map, authored explicitly so it can be argued with.

Ivan, 2026-08-07: 「整理好 EAR 三个都尽全力 compact ... 不损失结构性,中文都翻译成英语」.

**THE MEASUREMENT THAT JUSTIFIES THIS.** `#P16`'s own law says *an `A` that contains exactly one
`R` is a mis-cut arc*, because an arc is a decision made safe and a single belief update rarely
does that. Measured before touching anything:

    E01  10 arcs / 54 rounds   — every arc >= 2 rounds. WELL CUT, untouched.
    E02   5 arcs / 43 rounds   — one orphan (`A243`, numbered from a dead scheme).
    E03 117 arcs / 278 rounds  — **54 arcs hold exactly one round**. This is the mis-cut.

So the compaction is not a preference. E03 recorded 117 decisions; it made about 29.

**WHAT IS NOT LOST.** Round directories keep their identity and their `R` numbers — the law
requires `R` to be globally continuous and it already is. Nothing is deleted, nothing is renumbered
downward, no result moves. Only the ARC a round hangs under changes, and every arc name below is a
DECISION, stated in English.

**HOW THE GROUPS WERE CUT.** By the decision each thread made safe, read off the arc names
themselves — not by topic and not to hit a number. Where a thread genuinely made one decision over
many rounds it stays one arc however large; where fifteen arcs each recorded one step of a single
decision they become one.
"""

# ── E01 and E02: already well cut. Only the orphan moves. ────────────────────────────
KEEP_E01 = "every arc holds >= 2 rounds; nothing to merge"
E02_MERGE = {
    # `A243` is numbered from a scheme that no longer exists (`#911`: the sub-round scheme does not
    # exist) and holds one round. Its decision — are the walls I wrote measurable — is the same one
    # `A11`'s coupling audit closed.
    "A243_audit_the_walls": "A14_sexual_morality_versus_family_morality",
}

# ── E03: 117 -> 29. left = existing arc dir, right = the arc it becomes. ─────────────
E03_MERGE = {
    # ① IS SEXUAL MORALITY ONE OBJECT — the unity of the object, asked at three units.
    "A31_性到底是不是一件事":              "A31_is_sexual_morality_one_object",
    "A34_不靠池的那把尺子":                "A31_is_sexual_morality_one_object",
    "A37_社会层的性道德是不是一件事":       "A31_is_sexual_morality_one_object",
    "A38_人这个单位上性是不是也自成一条线":  "A31_is_sexual_morality_one_object",
    "A40_管别人的性是不是在管自己":         "A31_is_sexual_morality_one_object",

    # ② IS RESTRAINT THE ONLY OBJECT-GRADED AXIS — and why two instruments disagreed about it.
    "A32_性克制是不是唯一按对象分档的":      "A32_is_restraint_the_only_object_graded_axis",
    "A33_两具仪器为什么不一致":            "A32_is_restraint_the_only_object_graded_axis",
    "A36_那一族能不能机械化":              "A32_is_restraint_the_only_object_graded_axis",

    # ③ THE DECADE AS A UNIT, AND WHEN THE RIFT OPENED.
    "A35_年代这个单位":                    "A33_the_decade_as_a_unit_and_when_the_rift_opened",
    "A60_那五十年的扩大是不是匀速的":       "A33_the_decade_as_a_unit_and_when_the_rift_opened",
    "A61_九十年代那一格还没有解释":         "A33_the_decade_as_a_unit_and_when_the_rift_opened",
    "A64_裂开的时点是不是同一个":           "A33_the_decade_as_a_unit_and_when_the_rift_opened",
    "A66_裂还是漂要有一个连续的量":         "A33_the_decade_as_a_unit_and_when_the_rift_opened",
    "A67_把十年内全部年份用上":             "A33_the_decade_as_a_unit_and_when_the_rift_opened",
    "A68_放弃十年分辨率换回功效":           "A33_the_decade_as_a_unit_and_when_the_rift_opened",
    "A71_整张网格终于做多重性":             "A33_the_decade_as_a_unit_and_when_the_rift_opened",
    "A92_裂口自己在变宽是哪一边动的":       "A33_the_decade_as_a_unit_and_when_the_rift_opened",

    # ④ IS THE AXIS RELIGION OR POLITICS — the rival that would not die.
    "A39_那条弧会不会只是政治左右":         "A34_is_the_axis_religion_or_politics",
    "A65_宗教这个轴是不是可替换的":         "A34_is_the_axis_religion_or_politics",
    "A77_年龄和宗教是同一条缝吗":           "A34_is_the_axis_religion_or_politics",
    "A93_宗教是不是这个对象的那根轴":       "A34_is_the_axis_religion_or_politics",

    # ⑤ IS AN EXPLAINED SHARE A MEASUREMENT OR A DERIVATION.
    "A41_解释份额是测量还是推导":           "A35_is_an_explained_share_a_measurement_or_a_derivation",
    "A42_把页面上的份额改成推导":           "A35_is_an_explained_share_a_measurement_or_a_derivation",
    "A94_把两端拉开和解释了多少是同一件事吗": "A35_is_an_explained_share_a_measurement_or_a_derivation",
    "A97_一个差不能承载一个水平":           "A35_is_an_explained_share_a_measurement_or_a_derivation",

    # ⑥ DID THE DEVOUT CHANGE, AND IS IT A CEILING — the longest single decision in E03.
    "A43_两条已发表的话从没并排读过":       "A36_did_the_devout_change_or_is_it_a_ceiling",
    "A44_虔诚者改得少是不是天花板逼出来的":  "A36_did_the_devout_change_or_is_it_a_ceiling",
    "A45_虔诚者踩的是一个总闸还是一个个开关": "A36_did_the_devout_change_or_is_it_a_ceiling",
    "A47_最虔诚的三分之一是不是同一批人":    "A36_did_the_devout_change_or_is_it_a_ceiling",
    "A48_三堆是不是逐年取均值取出来的":      "A36_did_the_devout_change_or_is_it_a_ceiling",
    "A49_底堆是一堆人还是一片问不出来":      "A36_did_the_devout_change_or_is_it_a_ceiling",
    "A62_是人改了主意还是人换了":           "A36_did_the_devout_change_or_is_it_a_ceiling",
    "A63_虔诚者自己做了什么":              "A36_did_the_devout_change_or_is_it_a_ceiling",
    "A84_是世俗那边走了还是只在这一题上走了": "A36_did_the_devout_change_or_is_it_a_ceiling",
    "A89_虔诚有两张脸":                    "A36_did_the_devout_change_or_is_it_a_ceiling",

    # ⑦ HOW WIDE IS THE RIFT, AND DOES IT CLOSE.
    "A51_剩下那一对会不会也是小分母":       "A37_how_wide_is_the_rift_and_does_it_close",
    "A52_MDE关的是确认不是反驳":           "A37_how_wide_is_the_rift_and_does_it_close",
    "A53_顶对从哪儿到哪儿":                "A37_how_wide_is_the_rift_and_does_it_close",
    "A54_这条鸿沟在合还是在裂":            "A37_how_wide_is_the_rift_and_does_it_close",
    "A55_页面上最大的那个数一个对照都没有":  "A37_how_wide_is_the_rift_and_does_it_close",
    "A56_八条鸿沟每一条都欠同一个对照":     "A37_how_wide_is_the_rift_and_does_it_close",
    "A58_五裂三合一直只是方向":            "A37_how_wide_is_the_rift_and_does_it_close",
    "A59_八轮一具仪器":                    "A37_how_wide_is_the_rift_and_does_it_close",
    "A69_那个2点5倍是原始量表分":          "A37_how_wide_is_the_rift_and_does_it_close",
    "A70_那条幸存的事实自己稳不稳":         "A37_how_wide_is_the_rift_and_does_it_close",
    "A79_那一格里到底是谁在动":            "A37_how_wide_is_the_rift_and_does_it_close",
    "A80_缝到底有没有合上":                "A37_how_wide_is_the_rift_and_does_it_close",
    "A90_是错的和不许他教书是两件事吗":     "A37_how_wide_is_the_rift_and_does_it_close",
    "A95_那几块孤立的少数是不是同一批人":    "A37_how_wide_is_the_rift_and_does_it_close",
    "A96_不在任何一块里的反对者是谁":       "A37_how_wide_is_the_rift_and_does_it_close",

    # ⑧ AUDIT MY OWN WALLS, CORRECTIONS AND INSTRUMENTS — the self-audit thread, one decision:
    #    can the page be trusted to say what the evidence says.
    "A46_页面还按一个已知是错的切法组织着多少": "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A50_撤回要走到它引发的物件上":         "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A57_两条更正要走到页面上":            "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A72_两条更正落地":                    "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A73_那个匹配器的召回率":              "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A75_家族是什么这条轴从没扫过":         "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A76_那句话的主语从没被检验过":         "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A81_那三轮的总体到底虚高了多少":       "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A82_我写过的每一堵墙":                "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A83_墙倒了那条路走得通吗":            "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A85_那条线索经不经得住换一半人":       "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A86_不许说做不完先量一笔要多少":       "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A87_那条缝换一份问卷还在不在":         "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A88_我一轮前写下的那堵墙":            "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A98_还有几条守则只靠拷贝活着":         "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A99_那八轮的影响到底有多大":           "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A100_how_much_of_this_is_about_the_questionnaire": "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A101_what_can_still_be_rederived":     "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A102_does_the_evidence_still_say_what_the_ledger_says": "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A110_does_the_page_head_survive_its_own_audit": "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A113_the_repairs_that_lived_only_in_my_memory": "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A120_can_a_kill_be_made_to_fail":      "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A122_does_the_best_arc_survive_its_own_control": "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",
    "A133_the_register_is_a_search_result": "A38_can_the_page_be_trusted_to_say_what_the_evidence_says",

    # ⑨ THE SOCIETY AS THE UNIT.
    "A78_换单位到社会":                     "A39_the_society_as_the_unit",
    "A116_the_one_group_they_granted_less":  "A39_the_society_as_the_unit",
    "A117_the_society_unit_reopened":        "A39_the_society_as_the_unit",
    "A132_is_the_ordering_a_property_of_the_acts_or_of_the_society": "A39_the_society_as_the_unit",

    # ⑩ DOES THE COUPLING DEPEND ON THE ACT.
    "A103_does_the_coupling_depend_on_the_act":            "A40_does_the_coupling_depend_on_the_act",
    "A104_is_the_level_gap_about_acts_or_about_rulers":     "A40_does_the_coupling_depend_on_the_act",
    "A105_is_it_the_act_or_what_the_two_questions_ask":     "A40_does_the_coupling_depend_on_the_act",
    "A108_is_adultery_a_different_kind_of_judgement":       "A40_does_the_coupling_depend_on_the_act",

    # ⑪ WHICH INSTRUMENT CAN STILL BE ASKED.
    "A106_can_anything_i_hold_replicate_it":  "A41_which_instrument_can_still_be_asked",
    "A111_which_instrument_can_still_be_asked": "A41_which_instrument_can_still_be_asked",
    "A119_does_the_structure_exist_elsewhere": "A41_which_instrument_can_still_be_asked",

    # ⑫ ONE TIDE OR FOUR HISTORIES.
    "A107_one_tide_or_four_histories":                    "A42_one_tide_or_four_histories",
    "A109_is_not_one_tide_about_sex_or_about_surveys":     "A42_one_tide_or_four_histories",
    "A114_does_the_person_level_finding_explain_the_century": "A42_one_tide_or_four_histories",
    "A115_did_they_change_on_this_or_on_everything":       "A42_one_tide_or_four_histories",
    "A118_four_sexual_norms_one_scale":                    "A42_one_tide_or_four_histories",
    "A121_when_in_a_life_does_a_mind_change":              "A42_one_tide_or_four_histories",
    "A123_conversion_or_replacement_per_norm":             "A42_one_tide_or_four_histories",
    "A124_did_tolerance_grow_or_get_redirected":           "A42_one_tide_or_four_histories",

    # ⑬ IS THE PARTITION CONTENT OR MARGINALS — A125–A131, one decision, closed at `#959`.
    "A112_is_sexual_morality_a_distinct_object_at_the_person_level": "A43_is_the_partition_content_or_marginals",
    "A125_is_zero_sum_a_within_person_fact":               "A43_is_the_partition_content_or_marginals",
    "A126_is_the_two_layer_structure_general":             "A43_is_the_partition_content_or_marginals",
    "A127_is_the_partition_content_or_marginal":           "A43_is_the_partition_content_or_marginals",
    "A128_does_945_survive_a_marginal_null":               "A43_is_the_partition_content_or_marginals",
    "A129_does_the_residual_survive_a_full_marginal_null":  "A43_is_the_partition_content_or_marginals",
    "A130_the_partition_among_people_who_discriminate":     "A43_is_the_partition_content_or_marginals",
    "A131_why_the_contrast_reverses_among_sharp_discriminators": "A43_is_the_partition_content_or_marginals",
}

# A15–A30 are already cut at the decision level (2–9 rounds each) and keep their names.
KEEP_E03 = [f"A{n}" for n in range(15, 31)]
