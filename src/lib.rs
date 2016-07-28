extern crate rand;
use rand::Rng;

use std::cmp::min;

#[derive(Clone,Copy,PartialEq,Eq,Debug)]
pub struct Sample2D {
    x: usize,
    y: usize,
}

impl Sample2D {
    pub fn dstsqr(&self, other: &Sample2D) -> usize {
        let dx = (self.x as isize - other.x as isize).abs() as usize;
        let dy = (self.y as isize - other.y as isize).abs() as usize;
        dx*dx+dy*dy
    }
}

struct BackgroundGrid {
    data: Vec<usize>,
    dimensions: Vec<f64>,
    min_dst_sqr: f64,
    cell_size: usize,
    cell_count: Vec<usize>,
    cell_multiplicators: Vec<usize>,
}

impl BackgroundGrid {
    pub fn new(dimensions: Vec<f64>, min_distance: f64) -> BackgroundGrid {
        assert!(min_distance > 0.0);
        let dimension = dimensions.len();
        let cell_size = (min_distance / (dimension as f64).sqrt()).floor() as usize;
        let cell_count : Vec<usize> = dimensions.iter().map(|x| (x / (cell_size as f64)).ceil() as usize).collect();
        let data_size = cell_count.iter().fold(1_usize, |accu, x| accu * x);
        let mut cell_multiplicators = Vec::new();
        let mut multi_accu = 1_usize;
        for i in 0..dimension {
            cell_multiplicators.push(multi_accu);
            multi_accu *= cell_count[i];
        }
        BackgroundGrid {
            data: vec![0; data_size],
            dimensions: dimensions,
            min_dst_sqr: min_distance*min_distance,
            cell_size: cell_size,
            cell_count: cell_count,
            cell_multiplicators: cell_multiplicators,
        }
    }
    fn dst_sqr(x: &Vec<f64>, y: &Vec<f64>) -> f64 {
        debug_assert_eq!(x.len(), y.len());
        x.iter().zip(y.iter()).fold(1_f64, |accu, (xx, yx)| accu + xx*yx)
    }
    fn calc_idx(&self, cell_id: &Vec<usize>) -> usize {
        let mut idx=0;
        // TODO use multiplicators[0]=1 to optimize
        for i in 0..self.cell_multiplicators.len() {
            idx += cell_id[i] * self.cell_multiplicators[i];
        }
        idx
    }
    pub fn insert(&mut self, sample_position: Vec<f64>, samples: &mut Vec<Vec<f64>>) -> Result<usize,()> {
        if sample_position.iter().zip(self.dimensions.iter()).any(|(samp_x,dim_x)| samp_x >= dim_x) {
            return Err(());
        }
        let dimension = self.dimensions.len();
        assert_eq!(sample_position.len(), dimension);
        let cell_id : Vec<usize> = sample_position.iter().map(|x| *x as usize / self.cell_size).collect();
        let samp_idx = self.calc_idx(&cell_id);
        // TODO debug assert cell_id in range
        let cell_offs = (self.min_dst_sqr / (self.cell_size as f64)).ceil() as usize;
        let min_cell : Vec<usize> = cell_id.iter().map(|x| x.saturating_sub(cell_offs)).collect();
        let max_cell : Vec<usize> = cell_id.iter().zip(self.cell_count.iter()).map(|(x,size_x)| min(x + cell_offs, size_x-1)).collect();
        // TODO debug assert min_cell <= cell <= max_cell
        let mut indices = min_cell.clone();
        loop {
            let idx = self.calc_idx(&indices);
            match self.data[idx] {
                0 => (),
                other_id => {
                    let other_sample = &samples[other_id-1];
                    if BackgroundGrid::dst_sqr(&sample_position, other_sample) < self.min_dst_sqr {
                        return Err(());
                    }
                }
            }
            // iterate indices
            for i in 0..dimension {
                if indices[i] == max_cell[i] {
                    indices[i] = min_cell[i];
                } else {
                    indices[i]+=1;
                    break;
                }
            }
            // loop exit check
            if indices == max_cell {
                break;
            }
        }
        // no collission found
        samples.push(sample_position);
        assert_eq!(self.data[samp_idx], 0);
        self.data[samp_idx] = samples.len();
        Ok(samples.len())
    }
}

struct BackgroundGrid2D {
    data: Vec<Option<usize>>, // stores indices
    cells_x: usize, // for linearization of the 2D array
    cells_y: usize,
    width: usize,
    height: usize,
    min_distance_sqr: usize,
    cell_size: usize,
}

impl BackgroundGrid2D {
    pub fn new(width: usize, height: usize, min_distance: f64) -> BackgroundGrid2D {
        assert!(min_distance > 0.0_f64);
        let max_cell_size = (min_distance / 2.0_f64.sqrt()).floor();
        let cells_x : usize = (width as f64 / max_cell_size).ceil() as usize;
        let cells_y : usize = (height as f64 / max_cell_size).ceil() as usize;
        BackgroundGrid2D {
            data: vec![None; cells_x*cells_y],
            cells_x: cells_x,
            cells_y: cells_y,
            width: width,
            height: height,
            min_distance_sqr: (min_distance * min_distance).ceil() as usize,
            cell_size: max_cell_size as usize,
        }
    }
    pub fn insert(&mut self, sample_x: usize, sample_y: usize, samples: &mut Vec<Sample2D>) -> Result<usize,()> {
        let new_sample = Sample2D { x: sample_x, y: sample_y };
        if sample_x >= self.width || sample_y >= self.height {
            return Err(());
        }
        let cell_id_x = sample_x / self.cell_size;
        let cell_id_y = sample_y / self.cell_size;
        assert!(cell_id_x < self.cells_x);
        assert!(cell_id_y < self.cells_y);
        let min_cell_id_x = cell_id_x.saturating_sub(2);
        let max_cell_id_x = min((cell_id_x + 2),self.cells_x-1);
        let min_cell_id_y = cell_id_y.saturating_sub(2);
        let max_cell_id_y = min((cell_id_y + 2),self.cells_y-1);
        debug_assert!(min_cell_id_x <= cell_id_x);
        debug_assert!(max_cell_id_x >= cell_id_x);
        debug_assert!(min_cell_id_y <= cell_id_y);
        debug_assert!(max_cell_id_y >= cell_id_y);
        for y in min_cell_id_y..max_cell_id_y+1 {
            for x in min_cell_id_x..max_cell_id_x+1 {
                match self.data[y*self.cells_x+x] {
                    None => (),
                    Some(other_id) => {
                        let other_sample = &samples[other_id];
                        let dstsqr = new_sample.dstsqr(other_sample);
                        if dstsqr < self.min_distance_sqr {
                            return Err(());
                        }
                    }
                }
            }
        }
        // no collission found
        samples.push(new_sample);
        debug_assert_eq!(self.data[cell_id_y*self.cells_x+cell_id_x], None);
        self.data[cell_id_y*self.cells_x+cell_id_x] = Some(samples.len()-1);
        debug_assert_eq!(samples[samples.len()-1], new_sample);
        Ok(samples.len()-1)
    }
}

pub fn blue_noise (dimensions: Vec<f64>, min_distance: f64, k_abort: usize) -> Vec<Vec<f64>> {
    let mut rng = rand::thread_rng();
    let initial_sample = dimensions.iter().map(|x| rng.gen_range::<f64>(0_f64,*x)).collect();
    let mut samples : Vec<Vec<f64>> = Vec::new();
    let mut bggrid = BackgroundGrid::new(dimensions, min_distance);
    let initial_sample_id = bggrid.insert(initial_sample, &mut samples).unwrap();
    debug_assert_eq!(initial_sample_id, 1);
    let mut active : Vec<usize> = vec![initial_sample_id];
    while !active.is_empty() {
        let mut next_active : Vec<usize> = Vec::new();
        for current_id in active.iter() {
            let mut found_something = false;
            for _ in 0..k_abort {
                
            }
        }
    }
    samples
}

pub fn blue_noise2D (dimensions: [usize; 2], min_distance: f64, k_abort: usize) -> Vec<Sample2D> {
    let mut rng = rand::thread_rng();
    let initial_sample = Sample2D {
        x: rng.gen::<usize>() % dimensions[0],
        y: rng.gen::<usize>() % dimensions[1],
    };
    let mut samples : Vec<Sample2D> = Vec::new();
    let mut bggrid = BackgroundGrid2D::new(dimensions[0], dimensions[1], min_distance);
    let initial_sample_id = bggrid.insert(initial_sample.x, initial_sample.y, &mut samples).unwrap();
    debug_assert_eq!(initial_sample_id, 0);
    let mut active : Vec<usize> = vec![initial_sample_id];
    while !active.is_empty() {
        let mut next_active = Vec::new();
        for current_id in active.iter() {
            let current_sample = samples[*current_id].clone();
            let mut found_something = false;
            for _ in 0..k_abort {
                let distance = rng.gen_range::<f64>(min_distance, 2_f64 * min_distance);
                let angle = rng.gen_range::<f64>(0_f64, 2_f64 * std::f64::consts::PI);
                let xoffs = angle.cos() * distance;
                let yoffs = angle.sin() * distance;
                let test_samp_x = current_sample.x as f64 + xoffs;
                let test_samp_y = current_sample.y as f64 + yoffs;
                if test_samp_x < 0. || test_samp_y < 0. {
                    continue;
                }
                match bggrid.insert(test_samp_x.round() as usize, test_samp_y.round() as usize, &mut samples) {
                    Ok(new_id) => {
                        next_active.push(new_id);
                        found_something = true;
                    }
                    Err(_) => {
                        // wait for the next iteration
                    }
                }
            }
            if found_something {
                next_active.push(*current_id)
            }
        }
        active = next_active;
    }
    samples
}

#[test]
fn grid_corners() {
    let mut grid = BackgroundGrid2D::new(35,9,4.0);
    let mut samples = Vec::new();
    println!("cell size: {}", grid.cell_size);
    println!("grid dimensions: {} x {}", grid.cells_x, grid.cells_y);
    assert_eq!(grid.insert(0,9,&mut samples), Err(()));
    assert_eq!(grid.insert(0,0,&mut samples), Ok(0));
    assert_eq!(grid.insert(34,0,&mut samples), Ok(1));
    assert_eq!(grid.insert(0,8,&mut samples), Ok(2));
    assert_eq!(grid.insert(34,8,&mut samples), Ok(3));
}

#[cfg(test)]
mod tests {
    use super::*;
    fn get_image(radius: f64, size: usize) -> Vec<Vec<bool>> {
        let samples = blue_noise2D([size,size], radius, 30);
        let mut image = vec![vec![false; size]; size];
        for s in samples {
            image[s.y][s.x] = true;
        }
        image
    }
    fn output_ppm(image: &Vec<Vec<bool>>) {
        let size = image.len();
        assert_eq!(image[0].len(),size);
        println!("P1");
        println!("{} {}", size, size);
        for y in 0..size {
            for x in 0..size {
                print!("{} ", if image[y][x] { "1" } else { "0" });
            }
            println!("");
        }
    }
    #[test]
    fn output_sanity() {
        let size : isize = 128;
        let radius : isize = 8;
        let image = get_image(radius as f64, size as usize);
        for y in 0..size {
            for x in 0..size {
                if image[y as usize][x as usize] {
                    for dy in -radius..radius {
                        for dx in -radius..radius {
                            if dx == 0 && dy == 0 {
                                continue;
                            }
                            if -dy > y as isize || -dx > x as isize {
                                continue;
                            }
                            image.get((y+dy) as usize).map(|line| match line.get((x+dx) as usize) {
                                        Some(&true) => {
                                            assert!(dx*dx+dy*dy >= radius*radius);
                                        },
                                        _ => {}
                            });
                        }
                    }
                }
            }
        }
    }
    #[test]
    #[ignore]
    fn ppm() {
        let image = get_image(8.0, 128);
        output_ppm(&image);
    }
}
