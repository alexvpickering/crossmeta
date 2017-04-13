# Sources:
# https://www.jax.org/jax-mice-and-services/find-and-order-jax-mice/biospecimens/tissues-and-organs
# https://en.wikipedia.org/wiki/List_of_organs_of_the_human_body
# https://en.wikipedia.org/wiki/List_of_distinct_cell_types_in_the_adult_human_body


sources <- list(
	Cardiovascular = c(
		'heart',
			'cardiomyocytes',
			'sinoatrial node',
			'atrioventricular node',
			'Purkinje fibers',
			'ventricle',
			'atrium',
		'arteries',
		'aorta',
		'veins',
		'capillaries',
			'pericytes',
		'endothelial cells',
		'myocytes',
		'whole blood',
		'plasma',
		'serum'),




   	Lymphatic = c(
   		'lymph',
   		'lymphatic vessel',
   		'lymph nodes',
   		'bone marrow',
   			'bone marrow mesenchymal stem cells',
   		'thymus',
   		'spleen',
   		'gut-associated lymphoid tissue',
   		'tonsils'),




   	Integumentary = c(
   		'skin',
   			'keratinocytes',
   			'trichocytes',
   			'fibroblasts',
   			'myoepithelial cells',
   			'melanocytes',
   			'lacrimal glands',
   			'ceruminous glands',
   			'eccrine glands',
   			'apocrine sweat gland',
   			'glands of Moll',
   			'sebaceous glands',
   			'olfactory glands',

   		'fat',
   		'subcutaneous fat',
   		'visceral fat',
   			'white fat cells',
   				'adipose mesenchymal stem cells',
   			'brown fat cells',

   		'mammary glands'),



  	Musculoskeletal = c(
  		'bone',
  			'osteoclasts',
  			'osteoblasts',

  		'cartilage',
  			'chondrocytes',
  		'ligament',
  		'tendons',

  		'skeletal muscle',
  			'myosatellite cells'),



   	Reproductive = c(
   		'ovaries',
   			'oogonium',
   			'oocyte',
   		'fallopian tubes',
	 	'uterus',
	  	'vagina',
	  	"Bartholin's gland",
	  	'vulva',
        'clitoris',
        'placenta',
        'testes',
        	'spermatid',
        	'spermatocyte',
        	'spermatogonium',
        	'spermatozoon',
        'epididymis',
        'vas deferens',
        'seminal vesicles',
        'prostate',
        'bulbourethral gland',
        'penis',
        'urethral glands',
        'scrotum'),




   	Digestive = c(
   		'salivary glands',
   			"Von Ebner's glands",

   		'tongue',
   		'pharynx',
   		'oesophagus',

   		'stomach',
   			'gastric chief cells',
   			'gastric parietal cells',

        'small intestine',
     		'intestinal globlet cells',
        	'paneth cells',
        	'duodenal glands',

        'large intestine',
        'liver',
        	'hepatic stellate cells',
        	'hepatocyte',
        'gallbladder',
        'bile duct',
        'mesentery',

        'pancreas',
        	'centroacinar cells'),




   	Urinary = c(
   		'kidneys',
   			'podocytes',

   		'ureters',
   		'bladder',
   		'urethra'),





   	Respiratory = c(
   		'pharynx',
   		'larynx',

   		'lungs',
   			'trachea',
   			'bronchi',
   			'bronchioles',
   				'club cells',
   			'alveoli',
   				'type I alveolar cells',
   				'type II alveolar cells',
			'respiratory goblet cells',

  		'diaphragm'),




   	Endocrine = c(
   		'hypothalamus',
   			'magnocellular oxytocin cells',
			'magnocellular vasopressin cells',

   		'pituitary',
   			'somatotropes',
   			'prolactin cells',
   			'thyrotropes',
   			'gonadotropes',
   			'corticotropes',
   			'melanotropes',

   		'pineal',

   		'thyroid',
   			'follicular cells (thyroxine, triiodothryonine)',
   			'parafollicular cell (calcitonin)',

   		'parathyroid',
   			'parathyroid chief cell (parathyroid hormone)',
   			'oxyphil cell',

   		'adrenals',
   			'zona glomerulosa cells (aldosterone)',
   			'zona fasciculata cells (glucocorticoids)',
   			'zona reticularis cells (androgen precursors)',
   			'chromaffin cells (epinephrine, norepinephrine)',

   		'pancreatic islets',
   			'alpha cells (glucagon)',
        	'beta cells (insulin)',
        	'delta cells (somatostatin)',
        	'epsilon cells (ghrelin)',
       		'gamma cells (pancreatic polypeptide)',

   		'G cells (gastrin)',
		'I cells (cholecystokinin)',
        'S cells (secretin)',
        'K cells (glucose-dependent insulinotropic peptide)',
        'enterochromaffin cells (serotinin, histamine)',
        'N cells (neurotensin)',
		'L cells (glucagon-like peptide-1/2, pancreatic peptide YY)',


		"macula densa cells",
		"mesangial cells",
		'juxtaglomerular cells (renin)',

		'Leydig cells (testosterone)'),




   	Immune = c(
   		'peripheral blood mononuclear cells',
   		'hematopoietic stem cell',

   		'lymphocytes',
   		'common lymphoid progenitors',
   			'Natural Killer cells',
   			'T cells',
        		'CD4+ T cells',
        		'CD8+ T cells',
        		'regulatory T cells',
        	'B cells',
        	'plasma cells',

        'common myeloid progenitors',
	   		'monocytes',
	   			'macrophages',
	   			'dendritic cells',
	        'granulocytes',
	        	'basophils',
	        	'neutrophils',
	        	'eosinophils',
	        	'mast cell'),




   	Nervous = c(
   		'brain',
   			'head direction cells',

   		'cerebrum',
   			'olfactory bulb',
	   		'cerebral cortex',
	   			'corpus callosum',
	   			'Martinotti cells',
	   			'chandelier cells',
	   			'Cajal–Retzius cells',
	   			'spindle neuron',
	   			'grid cells',
	   			'speed cells',
	   			'Betz cells',
	   		'hippocampus',
	   			'dentate gyrus',
	   			'place cell',
	   			'boundary cells',
	        'basal ganglia',
	        	'pallidum',
		        'substantia nigra',
		        	'substantia nigra (pars reticulata)',
		        	'substantia nigra (pars compacta)',
		        'striatum',
		        	'medium spiny neurons',

        'diencephalon',
        	'thalamus',
        	'epithalamus',
        	'subthalamus',

        'brainstem',
        	'midbrain',
        	'pons',
        	'medulla oblongata',

        'cerebellum',
			'basket cells',
			'stellate cells',
			'golgi cells',
			'granule cells',
			'Lugaro cells',
			'unipolar brush cells',

        'spinal cord',
        	'spinal interneuron',
        	'Renshaw cells',
        'choroid plexus',

        'nerve',

    	'glia',
    		'oligodendrocytes',
    		'astrocyes',
    		'microglia',
    		'ependymal cells',
    			'tanycytes',

    		'schwann cells',
    		'satellite glial cells'),




   	Sensory = c(
   		'eyes',
   			'rod cells',
   			'cone cells',
   		'cornea',
   		'iris',
   		'ciliary body',
   		'lens',
   		'retina',
        'optic nerve',
        'pinna',
        'eardrum',
        'ossicles',
        'cochlea',
        	'organ of Corti',
        	'inner hair cell',
			'outer hair cell',


        'vestibule of the ear',
        'olfactory epithelium',
        'taste buds',
        'tactile epithelial cells',

        'glomus type I cell',
        'glomus type II cell')
   	)



setwd("~/Documents/Batcave/GEO/crossmeta")
load("~/Documents/Batcave/GEO/crossmeta/R/sysdata.rda")
devtools::use_data(gpl_bioc, homologene, token, sources, internal = TRUE, overwrite = TRUE)

