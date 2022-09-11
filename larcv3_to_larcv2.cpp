#include "larcv3_to_larcv2.h"

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



void larcv3_to_larcv2::convert(int n_events, int n_skip){

    int max_events = larcv3_manager.get_n_entries();

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
            larcv2_manager.set_id(larcv3_manager.event_id().run(),
            larcv3_manager.event_id().subrun(), 
            larcv3_manager.event_id().event());
        larcv2_manager.save_entry();
    }
    larcv2_manager.finalize();

}



void larcv3_to_larcv2::initialize(){

    // Initialize the input manager:
    for (auto in_file : input_files)
        larcv3_manager.add_in_file(in_file);
    larcv3_manager.initialize();

    larcv2_manager = larcv::IOManager(larcv::IOManager::kWRITE);
    larcv2_manager.set_out_file(output_file);
    larcv2_manager.initialize();
}



void larcv3_to_larcv2::add_in_file(const char * filename){
    std::string s(filename);
    input_files.push_back(filename);
}

void larcv3_to_larcv2::set_out_file(const char * filename){
    std::string s(filename);
    output_file = filename;
}



void larcv3_to_larcv2::convert_event(size_t entry){

    std::cout << "Converting entry " << entry << std::endl;
    // Read the entry for 2:
    larcv3_manager.read_entry(entry);

    for ( auto & producer : larcv3_manager.producer_list("particle")){
        convert_particle(producer);                
    }
    for ( auto & producer : larcv3_manager.producer_list("sparse2d")){
        convert_sparse2d(producer);
    }
    for ( auto & producer : larcv3_manager.producer_list("sparse3d")){
        convert_sparse3d(producer);
    }
    for ( auto & producer : larcv3_manager.producer_list("cluster2d")){
        convert_cluster2d(producer);
    }
    for ( auto & producer : larcv3_manager.producer_list("cluster3d")){
        convert_cluster3d(producer);
    }
    for ( auto & producer : larcv3_manager.producer_list("image2d")){
        convert_image2d(producer);
    }


}


// Conversion functions:
void larcv3_to_larcv2::convert_image2d(std::string producer){

    // Get the image2d from the input and output file:
    larcv::EventImage2D  * output_image_2d  = (larcv::EventImage2D * ) larcv2_manager.get_data("image2d", producer);
    std::shared_ptr<larcv3::EventTensor2D> input_image_2d = std::dynamic_pointer_cast<larcv3::EventTensor2D> (larcv3_manager.get_data("image2d", producer));
    for ( auto & image : input_image_2d->as_vector()){

        // Create a larcv3 meta for this:
        //larcv::ImageMeta  meta;

        larcv::ImageMeta meta=larcv::ImageMeta(image.meta().min(0), image.meta().min(1),
                 image.meta().max(0), image.meta().max(1),
                 image.meta().number_of_voxels(0), image.meta().number_of_voxels(1), 
                 static_cast<larcv::ProjectionID_t>(image.meta().id()),static_cast<larcv::DistanceUnit_t>(image.meta().unit()));

            // print()
        // meta.set_dimension(0, image.meta().width(), image.meta().rows(), image.meta().min()[0]);
        // meta.set_dimension(1, image.meta().height(), image.meta().cols(), image.meta().min()[1]);


        // Convert the input image to the output image by transposing:

        larcv::Image2D new_image(meta);

        float original_sum = 0;
        float new_sum = 0;
        // TODO: test this conversion!
        std::vector<size_t> coords;
        coords.resize(2);
        for (size_t i_row = 0; i_row < image.meta().rows(); i_row ++  ){
            for (size_t i_col = 0; i_col < image.meta().cols(); i_col ++ ){
                std::vector<size_t> rowcal{i_row, i_col};

                if (image.pixel(rowcal) != 0){
                    coords[0] = i_row;
                    coords[1] = i_col;
                    new_image.set_pixel(coords[0],coords[1], image.pixel(rowcal));
                    original_sum += image.pixel(rowcal);
                    new_sum += new_image.pixel(coords[0], coords[1]);
                }
            }
        }

        // Set the new image data:
        output_image_2d->emplace(std::move(new_image));
    }

}

void larcv3_to_larcv2::convert_particle(std::string producer){

    // Get the particles from the input and output file:
    larcv::EventParticle * output_particle = (larcv::EventParticle *) larcv2_manager.get_data("particle", producer);
    std::shared_ptr<larcv3::EventParticle> input_particle = std::dynamic_pointer_cast<larcv3::EventParticle>(larcv3_manager.get_data("particle", producer));

    for (auto & particle : input_particle->as_vector()){
        larcv::Particle new_particle;

        int shape = particle.shape();
        larcv::ShapeType_t new_shape = static_cast<larcv::ShapeType_t>(shape);

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
void larcv3_to_larcv2::convert_sparse2d(std::string producer){
    // std::cout << "Calling sparse2d for producer " << producer << std::endl;

    // Get the particles from the input and output file:
    larcv::EventSparseTensor2D  * output_sparse2d  = (larcv::EventSparseTensor2D *)  larcv2_manager.get_data("sparse2d", producer);
    std::shared_ptr<larcv3::EventSparseTensor2D> input_sparse2d = std::dynamic_pointer_cast<larcv3::EventSparseTensor2D> (larcv3_manager.get_data("sparse2d", producer));

    // print(producer, "Number of input sparse tensors: ", input_sparse2d.as_vector().size())
    for ( auto & sparse2d : input_sparse2d->as_vector()){
        
        // Create a larcv3 meta for this:
        //larcv::ImageMeta meta; 

        larcv::ImageMeta meta=larcv::ImageMeta(sparse2d.meta().min(0), sparse2d.meta().min(1),
                 sparse2d.meta().max(0), sparse2d.meta().max(1),
                 sparse2d.meta().number_of_voxels(0), sparse2d.meta().number_of_voxels(1), 
                 static_cast<larcv::ProjectionID_t>(sparse2d.meta().id()),static_cast<larcv::DistanceUnit_t>(sparse2d.meta().unit()));



        // meta.set_dimension(0, sparse2d.meta().width(),  sparse2d.meta().rows(), sparse2d.meta().min()[0]);
        // meta.set_dimension(1, sparse2d.meta().height(), sparse2d.meta().cols(), sparse2d.meta().min()[1]);
        // meta.set_projection_id(sparse2d.meta().id());

        // print(sparse2d.meta().dump())
        // Create a place to hold the output sparse2d:
        larcv::SparseTensor2D st;
        st.meta(meta);

        // Convert all of the voxels:
        for (auto & original_voxel : sparse2d.as_vector()){
            // std::vector<size_t> vec_of_coords;
            // vec_of_coords.push_back(sparse2d.meta().coordinates(original_voxel.id())[0]);
            // vec_of_coords.push_back(sparse2d.meta().coordinates(original_voxel.id())[1]);
            // auto new_index = meta.index(vec_of_coords[0],vec_of_coords[1]);
            // st.emplace(larcv::Voxel(new_index, original_voxel.value()));
            st.emplace(meta.pos_x(original_voxel.id()), meta.pos_y(original_voxel.id()),original_voxel.value());
                }

        // Set the new image data:
        output_sparse2d->emplace(std::move(st));
        // 
    }

    return;
}
void larcv3_to_larcv2::convert_sparse3d(std::string producer){
    // Get the tensors from the input and output file:
    larcv::EventSparseTensor3D  * output_sparse3d  = (larcv::EventSparseTensor3D * )  larcv2_manager.get_data("sparse3d", producer);
    std::shared_ptr<larcv3::EventSparseTensor3D> input_sparse3d = std::dynamic_pointer_cast<larcv3::EventSparseTensor3D> (larcv3_manager.get_data("sparse3d", producer));

    auto &original_meta = input_sparse3d->at(0).meta();

    // Create a larcv3 meta for this:
    // larcv::ImageMeta meta; 
    // meta.set_dimension(0, original_meta.width(),  original_meta.num_voxel_x(), original_meta.min[0]);
    // meta.set_dimension(1, original_meta.height(), original_meta.num_voxel_y(), original_meta.min[1]);
    // meta.set_dimension(2, original_meta.depth(), original_meta.num_voxel_z(), original_meta.min[2]);
    // meta.set_projection_id(0);

    larcv::Voxel3DMeta meta;

    meta.set(original_meta.min(0), original_meta.min(1), original_meta.min(2),
             original_meta.max(0), original_meta.max(1), original_meta.min(2),
             original_meta.number_of_voxels(0), original_meta.number_of_voxels(1), original_meta.number_of_voxels(2),
             static_cast<larcv::DistanceUnit_t>(original_meta.unit()));

        larcv::SparseTensor3D st;
    st.meta(meta);

    for (auto & original_voxel : input_sparse3d->as_vector()){
        // Convert all of the voxels:
        // std::vector<size_t> vec_of_coords;
        // vec_of_coords.push_back(original_meta.id_to_x_index(original_voxel.id()));
        // vec_of_coords.push_back(original_meta.id_to_y_index(original_voxel.id()));
        // vec_of_coords.push_back(original_meta.id_to_z_index(original_voxel.id()));
        // auto new_index = meta.index(vec_of_coords[0], vec_of_coords[1],vec_of_coords[2]);
        // st.emplace(larcv::Voxel(new_index, original_voxel.value()));
        st.emplace(meta.pos_x(original_voxel.id()), meta.pos_y(original_voxel.id()), meta.pos_z(original_voxel.id()), original_voxel.sum());
    }
        


    // Set the new image data:
    output_sparse3d->emplace(std::move(st),meta);

}
void larcv3_to_larcv2::convert_cluster2d(std::string producer){
    // Get the clusters from the input and output file:
larcv::EventClusterPixel2D  * output_cluster_2D  = (larcv::EventClusterPixel2D *)  larcv2_manager.get_data("cluster2d", producer);
    std::shared_ptr<larcv3::EventSparseCluster2D> input_cluster_2D = std::dynamic_pointer_cast<larcv3::EventSparseCluster2D> (larcv3_manager.get_data("cluster2d", producer));

    for (auto & cluster2d_set : input_cluster_2D->as_vector()){

        auto & original_meta = cluster2d_set.meta();

        // Create a larcv3 meta for this:
        // larcv::ImageMeta meta; 
        // meta.set_dimension(0, original_meta.width(),  original_meta.rows(), original_meta.min[0]);
        // meta.set_dimension(1, original_meta.height(), original_meta.cols(), original_meta.min[1]);
        // meta.set_projection_id(original_meta.id());


        larcv::ImageMeta meta=larcv::ImageMeta(original_meta.min(0), original_meta.min(1),
                 original_meta.max(0), original_meta.max(1),
                 original_meta.number_of_voxels(0), original_meta.number_of_voxels(1), 
                 static_cast<larcv::ProjectionID_t>(original_meta.id()),static_cast<larcv::DistanceUnit_t>(original_meta.unit()));

        // Create a place to hold the output cluster2d:
        larcv::ClusterPixel2D output_cluster2d_set;
        output_cluster2d_set.meta(meta);
        for (auto & cluster : cluster2d_set.as_vector()){
            //holder for new cluster:
            larcv::VoxelSet vs;
            vs.id(cluster.id());
            // Convert all of the voxels:
            for (auto & original_voxel : cluster.as_vector()){
                std::vector<size_t> vec_of_coords;
                vec_of_coords.push_back(original_meta.coordinates(original_voxel.id())[0]);
                vec_of_coords.push_back(original_meta.coordinates(original_voxel.id())[1]);
                auto new_index = meta.index(vec_of_coords[0],vec_of_coords[1]);
                vs.emplace(new_index, original_voxel.value(), false);
            }
            output_cluster2d_set.emplace(std::move(vs));
        }
        
           

        // Set the new image data:
        output_cluster_2D->emplace(std::move(output_cluster2d_set));
    }



}
void larcv3_to_larcv2::convert_cluster3d(std::string producer){
   
    larcv::EventClusterVoxel3D  * output_cluster_3D  = (larcv::EventClusterVoxel3D *)  larcv2_manager.get_data("cluster3d", producer);
    std::shared_ptr<larcv3::EventSparseCluster3D> input_cluster_3D = std::dynamic_pointer_cast<larcv3::EventSparseCluster3D> (larcv3_manager.get_data("cluster3d", producer));

    auto &original_meta = input_cluster_3D->at(0).meta();

    // Create a larcv3 meta for this:
    // larcv::ImageMeta meta; 
    // meta.set_dimension(0, original_meta.width(),  original_meta.num_voxel_x(), original_meta.min[0]);
    // meta.set_dimension(1, original_meta.height(), original_meta.num_voxel_y(), original_meta.min[1]);
    // meta.set_dimension(2, original_meta.depth(), original_meta.num_voxel_z(), original_meta.min[2]);
    // meta.set_projection_id(0);

    larcv::Voxel3DMeta meta;

    meta.set(original_meta.min(0), original_meta.min(1), original_meta.min(2),
             original_meta.max(0), original_meta.max(1), original_meta.min(2),
             original_meta.number_of_voxels(0), original_meta.number_of_voxels(1), original_meta.number_of_voxels(2),
             static_cast<larcv::DistanceUnit_t>(original_meta.unit()));

        // Create a place to hold the output cluster3d:
        larcv::ClusterVoxel3D output_cluster3d_set;
    output_cluster3d_set.meta(meta);


    for ( auto & cluster3d_set : input_cluster_3D->as_vector()){


        //holder for new cluster:
        larcv::VoxelSetArray vs;
        
        for (auto & original_voxel :  cluster3d_set.as_vector()){
            if (original_voxel.id() > meta.size()) continue;
            //vs.id(original_voxel.id());
            // Convert all of the voxels:
            // std::vector<size_t> vec_of_coords;
            // vec_of_coords.push_back(original_meta.coordinates(original_voxel.id())[0]);
            // vec_of_coords.push_back(original_meta.coordinates(original_voxel.id())[1]);
            // vec_of_coords.push_back(original_meta.coordinates(original_voxel.id())[2]);
            // auto new_index = meta.index(vec_of_coords[0],vec_of_coords[1],vec_of_coords[2]);
            larcv::VoxelSet vss;
            for (auto & avoxel :  original_voxel.as_vector())
            {
                vss.insert(larcv::Voxel(avoxel.id(), avoxel.value()));
            }

            vs.emplace(std::move(vss));
        }
        output_cluster3d_set.emplace(std::move(vs),meta);

    }

    // Set the new image data:
    output_cluster_3D->emplace(std::move(output_cluster3d_set),meta);
        
}
