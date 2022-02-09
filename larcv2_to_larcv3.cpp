#include "larcv2_to_larcv3.h"

// Larcv3 includes:
#include "larcv3/core/dataformat/EventTensor.h"
#include "larcv3/core/dataformat/EventParticle.h"
#include "larcv3/core/dataformat/EventSparseTensor.h"
#include "larcv3/core/dataformat/EventSparseCluster.h"


// Larcv2 includes:
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel2D.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"


void larcv2_to_larcv3::convert(int n_events, int n_skip){

    int max_events = larcv2_manager.get_n_entries();

    if (n_events == -1){
        n_events = max_events;
    }

    if (n_skip >= max_events) {
        std::cout << "Nothing to do, n_skip == " << n_skip 
                  << " < max_events == " << max_events << std::endl;
        return;
    }

    int range = n_events;
    // if there aren't enough events, reduce the range:
    if (range > max_events - n_skip){
        range = max_events - n_skip;
    }

    std::cout << "Range: " << range << std::endl;
    std::cout << "n_skip: " << n_skip << std::endl;

    
    for (size_t i = 0; i < range; i ++ ){
        convert_event(i+n_skip);
            larcv3_manager.set_id(larcv2_manager.event_id().run(),
            larcv2_manager.event_id().subrun(), 
            larcv2_manager.event_id().event());
        larcv3_manager.save_entry();
    }
    larcv3_manager.finalize();

}

void larcv2_to_larcv3::initialize(){

    // Initialize the input manager:
    for (auto in_file : input_files)
        larcv2_manager.add_in_file(in_file);
    larcv2_manager.initialize();

    larcv3_manager = larcv3::IOManager(larcv3::IOManager::kWRITE);
    larcv3_manager.set_out_file(output_file);
    larcv3_manager.initialize();
}

void larcv2_to_larcv3::add_in_file(const char * filename){
    std::string s(filename);
    input_files.push_back(filename);
}

void larcv2_to_larcv3::set_out_file(const char * filename){
    std::string s(filename);
    output_file = filename;
}


void larcv2_to_larcv3::convert_event(size_t entry){

    std::cout << "Converting entry " << entry << std::endl;
    // Read the entry for 2:
    larcv2_manager.read_entry(entry);

    for ( auto & producer : larcv2_manager.producer_list("particle")){
        convert_particle(producer);                
    }
    for ( auto & producer : larcv2_manager.producer_list("sparse2d")){
        convert_sparse2d(producer);
    }
    for ( auto & producer : larcv2_manager.producer_list("sparse3d")){
        convert_sparse3d(producer);
    }
    for ( auto & producer : larcv2_manager.producer_list("cluster2d")){
        convert_cluster2d(producer);
    }
    for ( auto & producer : larcv2_manager.producer_list("cluster3d")){
        convert_cluster3d(producer);
    }
    for ( auto & producer : larcv2_manager.producer_list("image2d")){
        convert_image2d(producer);
    }


}


// Conversion functions:
void larcv2_to_larcv3::convert_image2d(std::string producer){

    // Get the image2d from the input and output file:
    larcv::EventImage2D  * input_image_2d  = (larcv::EventImage2D * ) larcv2_manager.get_data("image2d", producer);
    larcv3::EventImage2D * output_image_2d = (larcv3::EventImage2D *) larcv3_manager.get_data("image2d", producer);

    for ( auto & image : input_image_2d->as_vector()){

        // Create a larcv3 meta for this:
        larcv3::ImageMeta2D  meta;
        
        // print()
        meta.set_dimension(0, image.meta().width(),  image.meta().rows(), image.meta().min_x());
        meta.set_dimension(1, image.meta().height(), image.meta().cols(), image.meta().min_y());
        meta.set_projection_id(image.meta().id());


        // Convert the input image to the output image by transposing:

        larcv3::Image2D new_image(meta);

        float original_sum = 0;
        float new_sum = 0;
        // TODO: test this conversion!
        std::vector<size_t> coords;
        coords.resize(2);
        for (size_t i_row = 0; i_row < image.meta().rows(); i_row ++  ){
            for (size_t i_col = 0; i_col < image.meta().cols(); i_col ++ ){
                if (image.pixel(i_row, i_col) != 0){
                    coords[0] = i_row;
                    coords[1] = i_col;
                    new_image.set_pixel(coords, image.pixel(i_row, i_col));
                    original_sum += image.pixel(i_row, i_col);
                    new_sum += new_image.pixel(coords);
                }
            }
        }

        // Set the new image data:
        output_image_2d->emplace(std::move(new_image));

    }

}

void larcv2_to_larcv3::convert_particle(std::string producer){

    // Get the particles from the input and output file:
    larcv::EventParticle * input_particle = (larcv::EventParticle *) larcv2_manager.get_data("particle", producer);
    larcv3::EventParticle * output_particle = (larcv3::EventParticle *) larcv3_manager.get_data("particle", producer);

    for (auto & particle : input_particle->as_vector()){
        larcv3::Particle new_particle;

        int shape = particle.shape();
        larcv3::ShapeType_t new_shape = static_cast<larcv3::ShapeType_t>(shape);

        new_particle.id(particle.id());
        new_particle.mcst_index(particle.mcst_index());
        new_particle.mct_index(particle.mct_index());
        new_particle.shape(new_shape);
        new_particle.nu_current_type(particle.nu_current_type());
        new_particle.nu_interaction_type(particle.nu_interaction_type());
        new_particle.track_id(particle.track_id());
        new_particle.pdg_code(particle.pdg_code());
        new_particle.momentum(particle.px(), particle.py(), particle.pz());
        new_particle.position(
            particle.position().x(), 
            particle.position().y(), 
            particle.position().z(), 
            particle.position().t());
        new_particle.end_position(
            particle.end_position().x(), 
            particle.end_position().y(), 
            particle.end_position().z(), 
            particle.end_position().t());
        new_particle.first_step(
            particle.first_step().x(), 
            particle.first_step().y(), 
            particle.first_step().z(), 
            particle.first_step().t());
        new_particle.last_step(
            particle.last_step().x(), 
            particle.last_step().y(), 
            particle.last_step().z(), 
            particle.last_step().t());
        new_particle.distance_travel(particle.distance_travel());
        new_particle.energy_init(particle.energy_init());
        new_particle.energy_deposit(particle.energy_deposit());
        new_particle.creation_process(particle.creation_process());
        new_particle.parent_track_id(particle.parent_track_id());
        new_particle.parent_pdg_code(particle.parent_pdg_code());
        new_particle.parent_position(
            particle.parent_position().x(), 
            particle.parent_position().y(), 
            particle.parent_position().z(), 
            particle.parent_position().t());
        new_particle.ancestor_track_id(particle.ancestor_track_id());
        new_particle.ancestor_pdg_code(particle.ancestor_pdg_code());
        new_particle.ancestor_position(
            particle.ancestor_position().x(), 
            particle.ancestor_position().y(), 
            particle.ancestor_position().z(), 
            particle.ancestor_position().t());

        // # BBox is not really supported yet
        // # boundingbox_2d
        // # boundingbox_3d

        output_particle->emplace_back(std::move(new_particle));
    }
}
void larcv2_to_larcv3::convert_sparse2d(std::string producer){
    // std::cout << "Calling sparse2d for producer " << producer << std::endl;

    // Get the particles from the input and output file:
    larcv::EventSparseTensor2D  * input_sparse2d  
        = (larcv::EventSparseTensor2D *)  larcv2_manager.get_data("sparse2d", producer);
    larcv3::EventSparseTensor2D * output_sparse2d 
        = (larcv3::EventSparseTensor2D * ) larcv3_manager.get_data("sparse2d", producer);

    // print(producer, "Number of input sparse tensors: ", input_sparse2d.as_vector().size())
    for ( auto & sparse2d : input_sparse2d->as_vector()){
        
        // Create a larcv3 meta for this:
        larcv3::ImageMeta2D meta; 
        meta.set_dimension(0, sparse2d.meta().width(),  sparse2d.meta().rows(), sparse2d.meta().min_x());
        meta.set_dimension(1, sparse2d.meta().height(), sparse2d.meta().cols(), sparse2d.meta().min_y());
        meta.set_projection_id(sparse2d.meta().id());

        // print(sparse2d.meta().dump())
        // Create a place to hold the output sparse2d:
        larcv3::SparseTensor2D st;
        st.meta(meta);

        // Convert all of the voxels:
        for (auto & original_voxel : sparse2d.as_vector()){
            std::vector<size_t> vec_of_coords;
            vec_of_coords.push_back(sparse2d.meta().index_to_row(original_voxel.id()));
            vec_of_coords.push_back(sparse2d.meta().index_to_col(original_voxel.id()));
            auto  new_index = meta.index(vec_of_coords);
            st.emplace(larcv3::Voxel(new_index, original_voxel.value()));
        }

        // Set the new image data:
        output_sparse2d->emplace(std::move(st));
        // 
    }

    return;
}
void larcv2_to_larcv3::convert_sparse3d(std::string producer){
    // Get the tensors from the input and output file:
    larcv::EventSparseTensor3D  * input_sparse3d  
        = (larcv::EventSparseTensor3D * )  larcv2_manager.get_data("sparse3d", producer);
    larcv3::EventSparseTensor3D * output_sparse3d 
        = (larcv3::EventSparseTensor3D * ) larcv3_manager.get_data("sparse3d", producer);

    auto & original_meta = input_sparse3d->meta();

    // Create a larcv3 meta for this:
    larcv3::ImageMeta3D meta; 
    meta.set_dimension(0, original_meta.width(),  original_meta.num_voxel_x(), original_meta.min_x());
    meta.set_dimension(1, original_meta.height(), original_meta.num_voxel_y(), original_meta.min_y());
    meta.set_dimension(2, original_meta.depth(), original_meta.num_voxel_z(), original_meta.min_z());
    meta.set_projection_id(0);
    larcv3::SparseTensor3D st;
    st.meta(meta);

    for (auto & original_voxel : input_sparse3d->as_vector()){
        // Convert all of the voxels:
        std::vector<size_t> vec_of_coords;
        vec_of_coords.push_back(original_meta.id_to_x_index(original_voxel.id()));
        vec_of_coords.push_back(original_meta.id_to_y_index(original_voxel.id()));
        vec_of_coords.push_back(original_meta.id_to_z_index(original_voxel.id()));
        auto new_index = meta.index(vec_of_coords);
        st.emplace(larcv3::Voxel(new_index, original_voxel.value()));
    }
        


    // Set the new image data:
    output_sparse3d->emplace(std::move(st));

}
void larcv2_to_larcv3::convert_cluster2d(std::string producer){
    // Get the clusters from the input and output file:
    larcv::EventClusterPixel2D  * input_cluster_2D  
        = (larcv::EventClusterPixel2D *)  larcv2_manager.get_data("cluster2d", producer);
    larcv3::EventSparseCluster2D * output_cluster_2D 
        = (larcv3::EventSparseCluster2D * ) larcv3_manager.get_data("cluster2d", producer);

    for (auto & cluster2d_set : input_cluster_2D->as_vector()){

        auto & original_meta = cluster2d_set.meta();

        // Create a larcv3 meta for this:
        larcv3::ImageMeta2D meta; 
        meta.set_dimension(0, original_meta.width(),  original_meta.rows(), original_meta.min_x());
        meta.set_dimension(1, original_meta.height(), original_meta.cols(), original_meta.min_y());
        meta.set_projection_id(original_meta.id());

        // Create a place to hold the output cluster2d:
        larcv3::SparseCluster2D output_cluster2d_set;
        output_cluster2d_set.meta(meta);
        for (auto & cluster : cluster2d_set.as_vector()){
            //holder for new cluster:
            larcv3::VoxelSet vs;
            vs.id(cluster.id());
            // Convert all of the voxels:
            for (auto & original_voxel : cluster.as_vector()){
                std::vector<size_t> vec_of_coords;
                vec_of_coords.push_back(original_meta.index_to_row(original_voxel.id()));
                vec_of_coords.push_back(original_meta.index_to_col(original_voxel.id()));
                auto new_index = meta.index(vec_of_coords);
                vs.emplace(new_index, original_voxel.value(), false);
            }
            output_cluster2d_set.emplace(std::move(vs));
        }
        
           

        // Set the new image data:
        output_cluster_2D->emplace(std::move(output_cluster2d_set));
    }



}
void larcv2_to_larcv3::convert_cluster3d(std::string producer){
   
    larcv::EventClusterVoxel3D  * input_cluster_3D  
        = (larcv::EventClusterVoxel3D *)  larcv2_manager.get_data("cluster3d", producer);
    larcv3::EventSparseCluster3D * output_cluster_3D 
        = (larcv3::EventSparseCluster3D * ) larcv3_manager.get_data("cluster3d", producer);
        
    auto & original_meta = input_cluster_3D->meta();

    // Create a larcv3 meta for this:
    larcv3::ImageMeta3D meta; 
    meta.set_dimension(0, original_meta.width(),  original_meta.num_voxel_x(), original_meta.min_x());
    meta.set_dimension(1, original_meta.height(), original_meta.num_voxel_y(), original_meta.min_y());
    meta.set_dimension(2, original_meta.depth(), original_meta.num_voxel_z(), original_meta.min_z());
    meta.set_projection_id(0);

    // Create a place to hold the output cluster3d:
    larcv3::SparseCluster3D output_cluster3d_set;
    output_cluster3d_set.meta(meta);


    for ( auto & cluster3d_set : input_cluster_3D->as_vector()){


        //holder for new cluster:
        larcv3::VoxelSet vs;
        vs.id(cluster3d_set.id());
        for (auto & original_voxel :  cluster3d_set.as_vector()){
            if (original_voxel.id() > meta.total_voxels()) continue;
            // Convert all of the voxels:
            std::vector<size_t> vec_of_coords;
            vec_of_coords.push_back(original_meta.id_to_x_index(original_voxel.id()));
            vec_of_coords.push_back(original_meta.id_to_y_index(original_voxel.id()));
            vec_of_coords.push_back(original_meta.id_to_z_index(original_voxel.id()));
            auto new_index = meta.index(vec_of_coords);
            vs.emplace(new_index, original_voxel.value(), false);
           
        }
        output_cluster3d_set.emplace(std::move(vs));

    }

    // Set the new image data:
    output_cluster_3D->emplace(std::move(output_cluster3d_set));
        
}